# QCManyBody Parallel Execution Module

This module provides parallel execution capabilities for QCManyBody calculations, enabling efficient computation of many-body expansion calculations across multiple CPU cores or compute nodes.

## Features

- **Multiple Execution Backends:**
  - Sequential (reference implementation, no parallelism)
  - Multiprocessing (single-node, multi-core parallelism)
  - MPI (multi-node distributed parallelism, coming soon)

- **Intelligent Task Management:**
  - Priority-based scheduling
  - Cost estimation
  - Error handling and retries
  - Progress tracking

- **Fault Tolerance:**
  - Checkpointing for long calculations
  - Result caching
  - Graceful error handling

## Quick Start

### Basic Usage

```python
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Enable parallel execution with default settings (auto-detect cores)
result = ManyBodyComputer.from_manybodyinput(
    input_model,
    parallel=True
)
```

### Custom Worker Count

```python
# Use 4 CPU cores
result = ManyBodyComputer.from_manybodyinput(
    input_model,
    parallel=True,
    n_workers=4
)
```

### Advanced Configuration

```python
from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

# Configure executor
config = ExecutorConfig(
    n_workers=8,
    timeout_per_task=3600,  # 1 hour per task
    max_retries=3,
    checkpoint_file="/tmp/checkpoint.json",
    cache_dir="/tmp/qcmb_cache",
    log_level="DEBUG"
)

# Create executor
executor = MultiprocessingExecutor(config)

# Use with context manager
with executor:
    result = ManyBodyComputer.from_manybodyinput(
        input_model,
        executor=executor
    )
```

## Architecture

### Core Classes

#### `BaseParallelExecutor`
Abstract base class defining the executor interface. All executors implement:
- `initialize()`: Set up resources
- `execute(tasks)`: Run tasks in parallel
- `shutdown()`: Clean up resources

#### `ParallelTask`
Dataclass encapsulating a single QC calculation:
- Task ID and labels
- Molecule and QC specification
- Priority and cost estimates
- Dependencies

#### `TaskResult`
Dataclass storing task execution results:
- Success/failure status
- QC results (AtomicResult)
- Error information
- Timing metadata

### Executors

#### `SequentialExecutor`
- Executes tasks one at a time
- Reference implementation for testing
- No overhead, guaranteed correctness
- Use for debugging or small systems

#### `MultiprocessingExecutor`
- Parallel execution using `multiprocessing.Pool`
- True parallelism (separate processes)
- Suitable for single-node, multi-core systems
- Automatic load balancing
- ~2-4x speedup typical on 4-8 cores

#### `MPIExecutor` (Coming Soon)
- Distributed execution using MPI
- Scales to 100+ nodes
- For HPC clusters
- Fault-tolerant master-worker architecture

## Performance Guidelines

### When to Use Parallel Execution

✅ **Good candidates:**
- Medium to large systems (5+ fragments)
- Multiple n-body levels
- Expensive methods (CCSD(T), high-level correlation)
- Multiple BSSE treatments

❌ **Poor candidates:**
- Small systems (<5 fragments)
- Fast methods (HF, DFT) with small basis sets
- Systems with very few tasks

### Optimal Worker Count

General guidelines:
- **Small systems** (5-10 fragments): 2-4 workers
- **Medium systems** (10-20 fragments): 4-8 workers
- **Large systems** (20+ fragments): 8-16 workers

The optimal number depends on:
- Available CPU cores
- Memory per task
- I/O requirements
- QC program threading (set QC program to 1 thread if using multiple workers)

### Expected Speedup

Typical speedup with `n` workers:
- **Ideal:** `n`x speedup
- **Good:** `0.7-0.9 * n`x speedup
- **Actual:** `0.5-0.8 * n`x speedup (accounting for overhead)

Factors affecting speedup:
- Task granularity (more tasks = better parallelism)
- Load balancing (similar task sizes = better)
- Overhead (process creation, IPC, serialization)

## Examples

### Example 1: Water Pentamer

```python
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Create water pentamer
mol = Molecule(
    symbols=["O", "H", "H"] * 5,
    geometry=[...],  # 15 atoms
    fragments=[[0,1,2], [3,4,5], [6,7,8], [9,10,11], [12,13,14]]
)

# Set up calculation
mbin = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {
            "bsse_type": ["cp"],
            "max_nbody": 3
        },
        "specification": {
            "mp2/cc-pvdz": {
                "program": "psi4",
                "model": {"method": "mp2", "basis": "cc-pvdz"}
            }
        }
    }
)

# Run in parallel
result = ManyBodyComputer.from_manybodyinput(
    mbin,
    parallel=True,
    n_workers=8
)

print(f"Interaction energy: {result.return_result:.6f} Eh")
```

### Example 2: Multi-level Calculation

```python
# Multi-level with different methods for different n-body levels
mbin = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {
            "bsse_type": ["cp", "nocp"],
            "levels": {
                1: "ccsd(t)/cc-pvtz",
                2: "mp2/cc-pvdz",
                3: "hf/cc-pvdz"
            }
        },
        "specification": {
            "ccsd(t)/cc-pvtz": {
                "program": "psi4",
                "model": {"method": "ccsd(t)", "basis": "cc-pvtz"}
            },
            "mp2/cc-pvdz": {
                "program": "psi4",
                "model": {"method": "mp2", "basis": "cc-pvdz"}
            },
            "hf/cc-pvdz": {
                "program": "psi4",
                "model": {"method": "hf", "basis": "cc-pvdz"}
            }
        }
    }
)

# Parallel execution automatically distributes all tasks
from qcmanybody.parallel import ExecutorConfig

config = ExecutorConfig(
    n_workers=16,
    timeout_per_task=7200,  # 2 hours for expensive CCSD(T)
    checkpoint_file="water_pentamer_checkpoint.json"
)

result = ManyBodyComputer.from_manybodyinput(
    mbin,
    parallel=True,
    executor_config=config
)
```

### Example 3: Custom Progress Tracking

```python
from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

# Progress callback
def progress_callback(task_id, completed, total):
    pct = 100 * completed / total
    print(f"Progress: {completed}/{total} ({pct:.1f}%) - Completed: {task_id}")

# Create executor
config = ExecutorConfig(n_workers=4)
executor = MultiprocessingExecutor(config)

# Manually execute with progress tracking
with executor:
    # Prepare tasks (normally done internally)
    from qcmanybody.core import ManyBodyCore

    core = ManyBodyCore(...)
    tasks = prepare_parallel_tasks(core)

    # Execute with callback
    results = executor.execute(tasks, progress_callback=progress_callback)
```

## Testing

Run parallel module tests:

```bash
# All tests
pytest qcmanybody/parallel/tests/

# Specific test file
pytest qcmanybody/parallel/tests/test_base.py

# With coverage
pytest qcmanybody/parallel/tests/ --cov=qcmanybody.parallel --cov-report=html
```

## Development Status

| Feature | Status | Notes |
|---------|--------|-------|
| Base Infrastructure | ✅ Complete | Abstract classes, data models |
| SequentialExecutor | ✅ Complete | Reference implementation |
| MultiprocessingExecutor | ✅ Complete | Single-node parallelism |
| Progress Tracking | ✅ Complete | Callbacks and logging |
| Error Handling | ✅ Complete | Retries and graceful failures |
| Task Scheduling | 🚧 In Progress | Priority-based scheduling |
| Checkpointing | 🚧 In Progress | Save/resume calculations |
| Result Caching | 🚧 In Progress | Cache task results |
| MPIExecutor | 📋 Planned | Multi-node parallelism |
| Dynamic Load Balancing | 📋 Planned | Work stealing |

## Contributing

See the main [development plan](../../PARALLEL_EXECUTION_PLAN.md) for:
- Architecture overview
- Development roadmap
- Coding standards
- Testing guidelines

## License

Same as QCManyBody main package.
