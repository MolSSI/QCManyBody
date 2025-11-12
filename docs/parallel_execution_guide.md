# Parallel Execution User Guide

## Overview

QCManyBody supports parallel execution of many-body expansion calculations, allowing you to leverage multiple CPU cores to significantly reduce computation time. This guide covers everything you need to know to use parallel execution effectively.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Concepts](#concepts)
3. [Usage Patterns](#usage-patterns)
4. [Configuration](#configuration)
5. [Performance Tuning](#performance-tuning)
6. [Troubleshooting](#troubleshooting)

---

## Quick Start

### Python API

```python
from qcmanybody import ParallelManyBodyComputer
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

# Option 1: Simple - auto-detect CPU cores
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True
)

# Option 2: Specify worker count
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True,
    n_workers=4
)

# Option 3: Full control with custom executor
config = ExecutorConfig(
    n_workers=4,
    timeout_per_task=1800.0,  # 30 minutes per task
    max_retries=2
)
executor = MultiprocessingExecutor(config)
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    executor=executor
)
```

### CLI

```bash
# Enable parallel execution with auto-detect
qcmanybody run input.json --parallel

# Specify number of workers
qcmanybody run input.json --parallel --n-workers 4

# Configure via input file
qcmanybody run input.json  # Uses execution settings from file
```

**Input file example:**
```json
{
  "execution": {
    "parallel": true,
    "n_workers": 4,
    "executor_type": "multiprocessing",
    "timeout_per_task": 1800.0,
    "max_retries": 2
  }
}
```

---

## Concepts

### Executors

Executors manage how tasks are executed:

- **SequentialExecutor**: Runs tasks one at a time (single-threaded)
  - Use for: Debugging, small calculations, limited resources
  - Pros: Simple, reliable, minimal overhead
  - Cons: No speedup from parallelization

- **MultiprocessingExecutor**: Runs tasks in parallel using multiple processes
  - Use for: Production runs, multi-core systems
  - Pros: True parallelism, scales with CPU cores
  - Cons: Process overhead, requires pickleable data

- **MPIExecutor** *(planned for Milestone 6)*: Distributed execution across nodes
  - Use for: HPC clusters, very large calculations
  - Pros: Scales to hundreds of cores
  - Cons: Requires MPI setup, more complex

### Task Parallelization

In a many-body expansion, each individual quantum chemistry calculation becomes a **task**:

```
3-body calculation with 4 fragments:
- 4 monomers (1-body)
- 6 dimers (2-body)
- 4 trimers (3-body)
= 14 independent tasks → can run in parallel
```

The parallel executor distributes these tasks across available CPU cores:

```
4 CPU cores, 14 tasks:
┌────────┬────────┬────────┬────────┐
│ Core 1 │ Core 2 │ Core 3 │ Core 4 │
├────────┼────────┼────────┼────────┤
│ Task 1 │ Task 2 │ Task 3 │ Task 4 │
│ Task 5 │ Task 6 │ Task 7 │ Task 8 │
│ Task 9 │Task 10 │Task 11 │Task 12 │
│Task 13 │Task 14 │   -    │   -    │
└────────┴────────┴────────┴────────┘
```

---

## Usage Patterns

### Pattern 1: Development & Testing

For development and debugging, use sequential execution:

```python
from qcmanybody import ManyBodyComputer

# Sequential (default) - easier to debug
result = ManyBodyComputer.from_manybodyinput(mb_input)
```

### Pattern 2: Production Runs

For production, enable parallel execution with auto-detect:

```python
from qcmanybody import ParallelManyBodyComputer

# Auto-detect CPU cores
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True
)
```

### Pattern 3: HPC/Batch Systems

For HPC systems, explicitly set worker count to match allocated resources:

```python
import os
from qcmanybody import ParallelManyBodyComputer

# Use SLURM_CPUS_PER_TASK or similar
n_workers = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))

result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True,
    n_workers=n_workers
)
```

### Pattern 4: Memory-Constrained Systems

For systems with limited memory, reduce worker count:

```python
# Use fewer workers to reduce memory usage
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True,
    n_workers=2  # Even on 8-core system
)
```

### Pattern 5: Long-Running Calculations

For calculations that might fail, use retries and timeouts:

```python
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

config = ExecutorConfig(
    n_workers=4,
    timeout_per_task=3600.0,  # 1 hour timeout per task
    max_retries=3,  # Retry failed tasks up to 3 times
    log_level="INFO"
)

executor = MultiprocessingExecutor(config)
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    executor=executor
)
```

---

## Configuration

### ExecutorConfig Options

```python
from qcmanybody.parallel import ExecutorConfig

config = ExecutorConfig(
    # Worker count
    n_workers=None,  # None = auto-detect CPU count

    # Timeout settings
    timeout_per_task=3600.0,  # Seconds, default = 1 hour

    # Retry settings
    max_retries=2,  # Number of retry attempts, default = 2

    # Checkpointing (future)
    checkpoint_interval=10,  # Save every N tasks
    checkpoint_file="checkpoint.json",  # Checkpoint file path

    # Caching (future)
    cache_dir=None,  # Directory for result cache

    # Logging
    log_level="INFO"  # DEBUG, INFO, WARNING, ERROR
)
```

### Choosing n_workers

**Rule of thumb:** Start with number of physical cores

```python
import os

# Get physical core count (not hyperthreads)
n_cores = os.cpu_count() // 2  # Divide by 2 if hyperthreading

# Conservative: leave 1 core for system
n_workers = max(1, n_cores - 1)
```

**Factors to consider:**

1. **CPU cores:** More workers = more parallelism (up to core count)
2. **Memory:** Each worker needs memory for QC calculation
3. **I/O:** Too many workers can bottleneck disk/network
4. **Task count:** No benefit if workers > tasks

**Examples:**

```python
# Laptop (4 cores, 8 GB RAM)
n_workers = 2

# Workstation (16 cores, 64 GB RAM)
n_workers = 8

# HPC node (64 cores, 256 GB RAM)
n_workers = 32
```

### Timeout Configuration

Set timeout based on expected task duration:

```python
# Small molecules (HF/STO-3G): 1-10 seconds
timeout_per_task = 60.0  # 1 minute buffer

# Medium molecules (MP2/cc-pVDZ): 1-10 minutes
timeout_per_task = 1800.0  # 30 minutes

# Large molecules (CCSD/cc-pVTZ): hours
timeout_per_task = 14400.0  # 4 hours
```

**Timeout behavior:**
- Task exceeding timeout is terminated
- Counts as task failure
- Can be retried if `max_retries > 0`

### Retry Configuration

```python
# Conservative (default): retry failures
max_retries = 2  # Try up to 3 times total

# Aggressive: no retries (fail fast)
max_retries = 0

# Patient: many retries (for flaky resources)
max_retries = 5
```

---

## Performance Tuning

### Expected Speedup

**Ideal speedup** = number of workers (e.g., 4 workers → 4x faster)

**Actual speedup** is usually 60-80% of ideal due to:
- Process overhead
- Task imbalance (some tasks take longer)
- System contention

**Example benchmarks:**

| System | Workers | Tasks | Sequential | Parallel | Speedup | Efficiency |
|--------|---------|-------|------------|----------|---------|------------|
| Laptop | 2 | 10 | 100s | 60s | 1.7x | 85% |
| Workstation | 4 | 20 | 200s | 60s | 3.3x | 83% |
| Workstation | 8 | 20 | 200s | 35s | 5.7x | 71% |

*Note: Efficiency = Speedup / Workers × 100%*

### Optimization Tips

#### 1. Match Workers to Task Count

```python
# BAD: More workers than tasks (wasted resources)
n_tasks = 5
n_workers = 16  # 11 workers idle

# GOOD: Workers ≤ tasks
n_workers = min(5, cpu_count)
```

#### 2. Balance Task Sizes

Tasks of similar duration lead to better load balancing:

```python
# Group tasks by computational cost
# (future feature: automatic scheduling)
```

#### 3. Monitor Resource Usage

```bash
# While calculation is running, monitor:
htop  # CPU usage per core
free -h  # Memory usage
iostat  # Disk I/O
```

#### 4. Profile Your Calculation

```python
import time
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

# Test different worker counts
for n in [1, 2, 4, 8]:
    config = ExecutorConfig(n_workers=n)
    executor = MultiprocessingExecutor(config)

    start = time.time()
    result = ParallelManyBodyComputer.from_manybodyinput(
        mb_input, executor=executor
    )
    elapsed = time.time() - start

    print(f"{n} workers: {elapsed:.1f}s (speedup: {elapsed_1worker/elapsed:.2f}x)")
```

### Memory Considerations

Each worker needs memory for:
- QC program (e.g., Psi4: 100-500 MB base)
- Molecule data (small)
- Basis set matrices (scales with molecule size)

**Example memory calculation:**

```
Small molecule (10 atoms, cc-pVDZ):
- Psi4 base: 200 MB
- Basis matrices: 50 MB
- Total per worker: ~300 MB
- 4 workers: ~1.2 GB

Medium molecule (30 atoms, cc-pVTZ):
- Psi4 base: 200 MB
- Basis matrices: 500 MB
- Total per worker: ~1 GB
- 4 workers: ~4 GB
```

**If you run out of memory:**
1. Reduce `n_workers`
2. Use smaller basis set
3. Request more memory (HPC)

---

## Troubleshooting

### Issue: No Speedup from Parallelization

**Symptoms:** Parallel takes same time or longer than sequential

**Causes:**
1. Too few tasks
2. Tasks too small (overhead dominates)
3. I/O bottleneck
4. Wrong executor type

**Solutions:**

```python
# Check task count
print(f"Total tasks: {len(tasks)}")
# Need at least 2x worker count for good speedup

# Check task duration
# If tasks < 1 second each, overhead may dominate

# Try different worker counts
for n in [1, 2, 4]:
    # ... benchmark each
```

### Issue: Out of Memory

**Symptoms:**
```
MemoryError
Killed (OOM)
```

**Solutions:**

```python
# Reduce worker count
config = ExecutorConfig(n_workers=2)  # Instead of 8

# Or use sequential executor
from qcmanybody import ManyBodyComputer
result = ManyBodyComputer.from_manybodyinput(mb_input)
```

### Issue: Hanging or Timeout

**Symptoms:** Calculation never completes, or tasks timeout

**Causes:**
1. Task actually too slow
2. QC program crashed
3. Network/filesystem issue

**Solutions:**

```python
# Increase timeout
config = ExecutorConfig(
    timeout_per_task=7200.0,  # 2 hours
    max_retries=1
)

# Enable debug logging
config = ExecutorConfig(log_level="DEBUG")

# Check logs for failed tasks
# Look for patterns in failures
```

### Issue: Tasks Fail Randomly

**Symptoms:** Some tasks succeed, some fail inconsistently

**Causes:**
1. Resource contention
2. Flaky QC program
3. Numerical instability

**Solutions:**

```python
# Increase retries
config = ExecutorConfig(max_retries=5)

# Reduce worker count (less contention)
config = ExecutorConfig(n_workers=2)

# Try sequential to isolate issue
from qcmanybody import ManyBodyComputer
result = ManyBodyComputer.from_manybodyinput(mb_input)
```

### Issue: "Can't pickle" Error

**Symptoms:**
```
PicklingError: Can't pickle <object>
```

**Cause:** Multiprocessing requires all data to be serializable

**Solution:** This should not occur with normal QCManyBody usage. If it does:
1. Report as a bug
2. Try sequential executor as workaround

### Getting Help

1. **Check logs:** Set `log_level="DEBUG"` for detailed output
2. **Minimal example:** Create smallest failing example
3. **Report issue:** https://github.com/MolSSI/QCManyBody/issues
   - Include: Python version, QCManyBody version, input file, error message

---

## Examples

### Example 1: Water Trimer (Simple)

```python
from qcmanybody import ParallelManyBodyComputer, ManyBodyInput
from qcelemental.models import Molecule

# Water trimer
mol = Molecule(
    symbols=["O", "H", "H"] * 3,
    geometry=...,  # coordinates
    fragments=[[0,1,2], [3,4,5], [6,7,8]]
)

mb_input = ManyBodyInput(
    molecule=mol,
    driver="energy",
    specification={
        "program": "psi4",
        "method": "mp2",
        "basis": "cc-pvdz"
    },
    keywords={
        "max_nbody": 3,
        "bsse_type": ["cp"]
    }
)

# Parallel execution
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True,
    n_workers=4
)

print(f"Total energy: {result.ret_energy}")
```

### Example 2: Large System with Custom Config

```python
from qcmanybody import ParallelManyBodyComputer
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

# Heavy calculation with long-running tasks
config = ExecutorConfig(
    n_workers=8,
    timeout_per_task=7200.0,  # 2 hours per task
    max_retries=3,
    log_level="INFO"
)

executor = MultiprocessingExecutor(config)

result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    executor=executor
)
```

### Example 3: HPC Batch Script

```bash
#!/bin/bash
#SBATCH --job-name=qcmb_parallel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Load modules
module load python/3.11
module load psi4/1.8

# Activate environment
source ~/.venvs/qcmanybody/bin/activate

# Run with parallel execution
qcmanybody run input.json \
    --parallel \
    --n-workers $SLURM_CPUS_PER_TASK \
    --output results.json
```

---

## Summary

**Key takeaways:**

1. ✅ **Enable parallel** with `parallel=True` or `--parallel`
2. ✅ **Auto-detect works** for most cases (omit `n_workers`)
3. ✅ **Expected speedup** is 60-80% of ideal
4. ✅ **Monitor resources** to optimize worker count
5. ✅ **Use timeouts** for long-running calculations
6. ✅ **Enable retries** for production robustness

**Next steps:**
- Try parallel execution on your calculations
- Benchmark to find optimal worker count
- Report any issues or performance concerns

For more information:
- API Reference: `qcmanybody/parallel/README.md`
- Development Plan: `PARALLEL_EXECUTION_PLAN.md`
- Examples: `examples/cli/05_parallel_basic.json`
