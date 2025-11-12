# Migration Guide: Adding Parallel Execution to Existing Code

This guide helps existing QCManyBody users adopt parallel execution with minimal code changes.

## Table of Contents

1. [Overview](#overview)
2. [Backward Compatibility](#backward-compatibility)
3. [Migration Strategies](#migration-strategies)
4. [Code Examples](#code-examples)
5. [Common Patterns](#common-patterns)
6. [Testing Your Migration](#testing-your-migration)
7. [Troubleshooting](#troubleshooting)

---

## Overview

Parallel execution in QCManyBody is:

- ✅ **Backward compatible**: Existing code continues to work unchanged
- ✅ **Opt-in**: Enable only where you need it
- ✅ **Drop-in replacement**: Minimal code changes required
- ✅ **Gradual adoption**: Migrate incrementally, test thoroughly

**Philosophy:** "Make sequential easy, parallel straightforward, advanced possible"

---

## Backward Compatibility

### Guarantee

**All existing code continues to work unchanged.**

```python
# This still works exactly as before
from qcmanybody import ManyBodyComputer

result = ManyBodyComputer.from_manybodyinput(mb_input)
# Sequential execution (default behavior)
```

### No Breaking Changes

- Default behavior unchanged (sequential)
- All existing APIs preserved
- No deprecated features
- CLI remains compatible

---

## Migration Strategies

### Strategy 1: Quick Win (5 minutes)

**Goal:** Enable parallel execution with single line change

**Before:**

```python
from qcmanybody import ManyBodyComputer

result = ManyBodyComputer.from_manybodyinput(mb_input)
```

**After:**

```python
from qcmanybody import ParallelManyBodyComputer

result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True  # ← Only change needed
)
```

**Benefits:**
- ✅ Immediate speedup (typically 2-4x on 4-core systems)
- ✅ Auto-detects CPU cores
- ✅ Minimal risk

**When to use:**
- Development/testing on multi-core workstation
- Quick performance boost needed
- No specific tuning requirements

---

### Strategy 2: Controlled Migration (30 minutes)

**Goal:** Enable parallel with explicit configuration

**Before:**

```python
from qcmanybody import ManyBodyComputer

result = ManyBodyComputer.from_manybodyinput(mb_input)
```

**After:**

```python
from qcmanybody import ParallelManyBodyComputer
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

# Explicit configuration
config = ExecutorConfig(
    n_workers=4,
    timeout_per_task=1800.0,
    max_retries=2,
    log_level="INFO"
)

executor = MultiprocessingExecutor(config)

result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    executor=executor
)
```

**Benefits:**
- ✅ Full control over parallelization
- ✅ Explicit resource usage
- ✅ Production-ready configuration

**When to use:**
- Production environments
- HPC clusters (explicit core counts)
- Resource-constrained systems
- Long-running calculations

---

### Strategy 3: Feature Flag Migration (1-2 hours)

**Goal:** Add parallel execution as optional feature

**Before:**

```python
def run_calculation(mb_input):
    from qcmanybody import ManyBodyComputer
    result = ManyBodyComputer.from_manybodyinput(mb_input)
    return result
```

**After:**

```python
def run_calculation(mb_input, parallel=False, n_workers=None):
    if parallel:
        from qcmanybody import ParallelManyBodyComputer
        result = ParallelManyBodyComputer.from_manybodyinput(
            mb_input,
            parallel=True,
            n_workers=n_workers
        )
    else:
        from qcmanybody import ManyBodyComputer
        result = ManyBodyComputer.from_manybodyinput(mb_input)

    return result

# Usage
result = run_calculation(mb_input, parallel=True, n_workers=4)
```

**Benefits:**
- ✅ Gradual rollout
- ✅ Easy A/B testing
- ✅ Safe fallback to sequential

**When to use:**
- Large codebase with many call sites
- Gradual feature rollout
- Need for easy rollback

---

## Code Examples

### Example 1: Simple Script

**Before:**

```python
#!/usr/bin/env python
"""Simple many-body calculation."""

from qcmanybody import ManyBodyComputer
from qcelemental.models import Molecule

# Setup molecule
mol = Molecule(
    symbols=["He", "He", "He"],
    geometry=[[0, 0, 0], [0, 0, 3], [0, 3, 0]],
    fragments=[[0], [1], [2]]
)

# Create input
mb_input = create_manybody_input(mol, method="mp2", basis="cc-pvdz")

# Run calculation
result = ManyBodyComputer.from_manybodyinput(mb_input)

# Print results
print(f"Total energy: {result.ret_energy}")
```

**After (minimal change):**

```python
#!/usr/bin/env python
"""Simple many-body calculation with parallel execution."""

from qcmanybody import ParallelManyBodyComputer  # ← Changed import
from qcelemental.models import Molecule

# Setup molecule (unchanged)
mol = Molecule(
    symbols=["He", "He", "He"],
    geometry=[[0, 0, 0], [0, 0, 3], [0, 3, 0]],
    fragments=[[0], [1], [2]]
)

# Create input (unchanged)
mb_input = create_manybody_input(mol, method="mp2", basis="cc-pvdz")

# Run calculation (parallel)
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True  # ← Added parameter
)

# Print results (unchanged)
print(f"Total energy: {result.ret_energy}")
```

---

### Example 2: Function with Optional Parallel

**Before:**

```python
def compute_binding_energy(fragments, method="mp2", basis="cc-pvdz"):
    """Compute binding energy for given fragments."""
    from qcmanybody import ManyBodyComputer

    # Create molecule
    mol = create_molecule_from_fragments(fragments)

    # Create input
    mb_input = create_manybody_input(mol, method=method, basis=basis)

    # Compute
    result = ManyBodyComputer.from_manybodyinput(mb_input)

    return result.ret_energy
```

**After:**

```python
def compute_binding_energy(
    fragments,
    method="mp2",
    basis="cc-pvdz",
    parallel=False,  # ← New parameter
    n_workers=None   # ← New parameter
):
    """Compute binding energy for given fragments.

    Parameters
    ----------
    fragments : list
        Fragment geometries
    method : str
        QC method
    basis : str
        Basis set
    parallel : bool, optional
        Enable parallel execution
    n_workers : int, optional
        Number of workers (None = auto-detect)
    """
    # Select computer class based on parallel flag
    if parallel:
        from qcmanybody import ParallelManyBodyComputer as Computer
        kwargs = {"parallel": True, "n_workers": n_workers}
    else:
        from qcmanybody import ManyBodyComputer as Computer
        kwargs = {}

    # Create molecule (unchanged)
    mol = create_molecule_from_fragments(fragments)

    # Create input (unchanged)
    mb_input = create_manybody_input(mol, method=method, basis=basis)

    # Compute with appropriate executor
    result = Computer.from_manybodyinput(mb_input, **kwargs)

    return result.ret_energy

# Usage examples:
# Sequential (default, backward compatible)
energy = compute_binding_energy(fragments)

# Parallel with auto-detect
energy = compute_binding_energy(fragments, parallel=True)

# Parallel with explicit workers
energy = compute_binding_energy(fragments, parallel=True, n_workers=4)
```

---

### Example 3: Class-Based Application

**Before:**

```python
class ManyBodyCalculator:
    """Manage many-body calculations."""

    def __init__(self, method="mp2", basis="cc-pvdz"):
        self.method = method
        self.basis = basis

    def compute(self, molecule):
        """Run calculation for molecule."""
        from qcmanybody import ManyBodyComputer

        mb_input = self._create_input(molecule)
        result = ManyBodyComputer.from_manybodyinput(mb_input)

        return result

    def _create_input(self, molecule):
        """Create ManyBodyInput."""
        # ... implementation
```

**After:**

```python
class ManyBodyCalculator:
    """Manage many-body calculations with optional parallel execution."""

    def __init__(
        self,
        method="mp2",
        basis="cc-pvdz",
        parallel=False,      # ← New parameter
        n_workers=None,      # ← New parameter
        executor_config=None # ← New parameter
    ):
        self.method = method
        self.basis = basis
        self.parallel = parallel
        self.n_workers = n_workers
        self.executor_config = executor_config
        self._executor = None

    def compute(self, molecule):
        """Run calculation for molecule.

        Uses parallel execution if enabled in constructor.
        """
        mb_input = self._create_input(molecule)

        if self.parallel:
            from qcmanybody import ParallelManyBodyComputer

            # Use custom executor if provided
            if self.executor_config:
                from qcmanybody.parallel import MultiprocessingExecutor
                executor = MultiprocessingExecutor(self.executor_config)
                result = ParallelManyBodyComputer.from_manybodyinput(
                    mb_input,
                    executor=executor
                )
            else:
                result = ParallelManyBodyComputer.from_manybodyinput(
                    mb_input,
                    parallel=True,
                    n_workers=self.n_workers
                )
        else:
            from qcmanybody import ManyBodyComputer
            result = ManyBodyComputer.from_manybodyinput(mb_input)

        return result

    def _create_input(self, molecule):
        """Create ManyBodyInput (unchanged)."""
        # ... implementation

# Usage examples:

# Sequential (backward compatible)
calc = ManyBodyCalculator()
result = calc.compute(molecule)

# Parallel with auto-detect
calc = ManyBodyCalculator(parallel=True)
result = calc.compute(molecule)

# Parallel with explicit configuration
from qcmanybody.parallel import ExecutorConfig

config = ExecutorConfig(n_workers=4, timeout_per_task=1800.0)
calc = ManyBodyCalculator(parallel=True, executor_config=config)
result = calc.compute(molecule)
```

---

### Example 4: CLI Application

**Before:**

```python
def main(args):
    """Main CLI entry point."""
    from qcmanybody import ManyBodyComputer

    # Parse input
    mb_input = parse_input_file(args.input)

    # Run calculation
    result = ManyBodyComputer.from_manybodyinput(mb_input)

    # Write output
    write_output(result, args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("-o", "--output")

    args = parser.parse_args()
    main(args)
```

**After:**

```python
def main(args):
    """Main CLI entry point with parallel support."""
    # Parse input
    mb_input = parse_input_file(args.input)

    # Select computer based on args
    if args.parallel:
        from qcmanybody import ParallelManyBodyComputer
        result = ParallelManyBodyComputer.from_manybodyinput(
            mb_input,
            parallel=True,
            n_workers=args.n_workers
        )
    else:
        from qcmanybody import ManyBodyComputer
        result = ManyBodyComputer.from_manybodyinput(mb_input)

    # Write output
    write_output(result, args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("-o", "--output")

    # Add parallel arguments
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Enable parallel execution"
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=None,
        help="Number of workers (default: auto-detect)"
    )

    args = parser.parse_args()
    main(args)
```

---

## Common Patterns

### Pattern 1: Environment-Based Configuration

Enable parallel based on environment variable:

```python
import os
from qcmanybody import ManyBodyComputer, ParallelManyBodyComputer

# Check environment
parallel = os.environ.get("QCMANYBODY_PARALLEL", "false").lower() == "true"
n_workers = int(os.environ.get("QCMANYBODY_WORKERS", "0")) or None

# Select computer
if parallel:
    result = ParallelManyBodyComputer.from_manybodyinput(
        mb_input,
        parallel=True,
        n_workers=n_workers
    )
else:
    result = ManyBodyComputer.from_manybodyinput(mb_input)
```

**Usage:**

```bash
# Sequential
python script.py

# Parallel with auto-detect
QCMANYBODY_PARALLEL=true python script.py

# Parallel with 4 workers
QCMANYBODY_PARALLEL=true QCMANYBODY_WORKERS=4 python script.py
```

---

### Pattern 2: Configuration File

Read parallel settings from config file:

```python
import json
from qcmanybody import ManyBodyComputer, ParallelManyBodyComputer
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

def load_config(config_file):
    """Load configuration from JSON file."""
    with open(config_file) as f:
        return json.load(f)

def run_with_config(mb_input, config):
    """Run calculation with configuration."""
    parallel_config = config.get("parallel", {})

    if parallel_config.get("enabled", False):
        # Create executor config
        exec_config = ExecutorConfig(
            n_workers=parallel_config.get("n_workers"),
            timeout_per_task=parallel_config.get("timeout", 3600.0),
            max_retries=parallel_config.get("max_retries", 2),
            log_level=parallel_config.get("log_level", "INFO")
        )

        executor = MultiprocessingExecutor(exec_config)

        result = ParallelManyBodyComputer.from_manybodyinput(
            mb_input,
            executor=executor
        )
    else:
        result = ManyBodyComputer.from_manybodyinput(mb_input)

    return result

# Usage
config = load_config("config.json")
result = run_with_config(mb_input, config)
```

**config.json:**

```json
{
  "parallel": {
    "enabled": true,
    "n_workers": 4,
    "timeout": 1800.0,
    "max_retries": 2,
    "log_level": "INFO"
  }
}
```

---

### Pattern 3: Wrapper Function

Create convenience wrapper:

```python
def compute_manybody(
    mb_input,
    parallel="auto",
    n_workers=None,
    **kwargs
):
    """Compute many-body expansion with smart parallel selection.

    Parameters
    ----------
    mb_input : ManyBodyInput
        Input specification
    parallel : bool or "auto"
        - True: Always use parallel
        - False: Always use sequential
        - "auto": Use parallel if > 4 cores available
    n_workers : int, optional
        Number of workers
    **kwargs
        Additional arguments for ParallelManyBodyComputer

    Returns
    -------
    ManyBodyResult
    """
    # Auto-detection logic
    if parallel == "auto":
        import os
        parallel = os.cpu_count() > 4

    # Select and execute
    if parallel:
        from qcmanybody import ParallelManyBodyComputer
        result = ParallelManyBodyComputer.from_manybodyinput(
            mb_input,
            parallel=True,
            n_workers=n_workers,
            **kwargs
        )
    else:
        from qcmanybody import ManyBodyComputer
        result = ManyBodyComputer.from_manybodyinput(mb_input, **kwargs)

    return result

# Usage examples:
result = compute_manybody(mb_input)  # Auto-detect
result = compute_manybody(mb_input, parallel=True)  # Force parallel
result = compute_manybody(mb_input, parallel=False)  # Force sequential
```

---

## Testing Your Migration

### Step 1: Verify Correctness

Ensure parallel gives identical results to sequential:

```python
from qcmanybody import ManyBodyComputer, ParallelManyBodyComputer
import numpy as np

# Run sequential
result_seq = ManyBodyComputer.from_manybodyinput(mb_input)

# Run parallel
result_par = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True
)

# Compare results
assert np.isclose(result_seq.ret_energy, result_par.ret_energy, rtol=1e-10)
print("✓ Results match!")
```

### Step 2: Benchmark Performance

Measure actual speedup:

```python
import time

# Sequential
start = time.time()
result_seq = ManyBodyComputer.from_manybodyinput(mb_input)
seq_time = time.time() - start

# Parallel
start = time.time()
result_par = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True,
    n_workers=4
)
par_time = time.time() - start

# Report
speedup = seq_time / par_time
print(f"Sequential: {seq_time:.1f}s")
print(f"Parallel (4 workers): {par_time:.1f}s")
print(f"Speedup: {speedup:.2f}x")
```

### Step 3: Test Error Handling

Verify errors are caught gracefully:

```python
# Test with invalid configuration
try:
    result = ParallelManyBodyComputer.from_manybodyinput(
        mb_input,
        parallel=True,
        n_workers=0  # Invalid!
    )
except ValueError as e:
    print(f"✓ Caught invalid config: {e}")

# Test with timeout
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

config = ExecutorConfig(
    n_workers=2,
    timeout_per_task=1.0,  # Very short timeout
    max_retries=0
)

executor = MultiprocessingExecutor(config)
result = ParallelManyBodyComputer.from_manybodyinput(
    mb_input,
    executor=executor
)

# Check for timeouts
# (some tasks may timeout, check result status)
```

---

## Troubleshooting

### Issue: No Speedup After Migration

**Checklist:**

1. ✅ Verify parallel execution is actually enabled
   ```python
   # Add logging to confirm
   import logging
   logging.basicConfig(level=logging.INFO)
   ```

2. ✅ Check CPU usage during run
   ```bash
   htop  # Should see multiple cores at 100%
   ```

3. ✅ Ensure enough tasks for parallelization
   ```python
   # Need at least 2x worker count for good speedup
   print(f"Tasks: {count_tasks(mb_input)}")
   print(f"Workers: {n_workers}")
   ```

4. ✅ Rule out I/O bottleneck
   ```python
   # Try fewer workers
   result = ParallelManyBodyComputer.from_manybodyinput(
       mb_input,
       parallel=True,
       n_workers=2  # Instead of 8
   )
   ```

### Issue: Results Don't Match

**Possible causes:**

1. Numerical precision differences (normal, acceptable if small)
2. Race condition in QC program (rare)
3. Bug in parallel implementation (report it!)

**Debug:**

```python
# Run multiple times and compare
results = []
for i in range(5):
    result = ParallelManyBodyComputer.from_manybodyinput(
        mb_input,
        parallel=True
    )
    results.append(result.ret_energy)

# Check variance
import numpy as np
print(f"Mean: {np.mean(results)}")
print(f"Std: {np.std(results)}")  # Should be < 1e-9
```

### Issue: Migration Breaks Existing Tests

**Common causes:**

1. Tests check exact execution order
2. Tests mock `ManyBodyComputer` specifically
3. Tests validate log output format

**Solutions:**

```python
# Make tests flexible about computer class
from qcmanybody import ManyBodyComputer, ParallelManyBodyComputer

# Instead of:
# computer = ManyBodyComputer.from_manybodyinput(...)

# Use:
computer_class = get_computer_class()  # From config/fixture
computer = computer_class.from_manybodyinput(...)
```

---

## Best Practices

### ✅ Do:

- Start with small test case
- Verify results match before scaling up
- Add parallel as optional feature first
- Document new configuration options
- Test both sequential and parallel paths
- Monitor resource usage in production

### ❌ Don't:

- Force parallel for all use cases
- Remove sequential fallback
- Ignore backward compatibility
- Deploy without testing
- Use parallel for tiny calculations (overhead dominates)

---

## Migration Checklist

Use this checklist when migrating:

- [ ] Identify code locations that call `ManyBodyComputer`
- [ ] Determine which need parallel execution
- [ ] Add parallel parameters to function signatures
- [ ] Update imports where needed
- [ ] Add configuration options (CLI args, config file, etc.)
- [ ] Write tests for both sequential and parallel paths
- [ ] Benchmark to verify speedup
- [ ] Document new parallel options for users
- [ ] Update examples and tutorials
- [ ] Test in staging/development environment
- [ ] Deploy incrementally with monitoring
- [ ] Update production documentation

---

## Getting Help

**Questions?**

- Check: `docs/parallel_execution_guide.md`
- API reference: `docs/parallel_api_reference.md`
- Examples: `examples/cli/05_parallel_basic.json`

**Issues?**

- Report bugs: https://github.com/MolSSI/QCManyBody/issues
- Include: Minimal failing example, error message, system info

**Need features?**

- Feature requests welcome!
- Describe use case and requirements
- Contribute via pull request

---

## Summary

**Key points:**

1. ✅ Parallel execution is fully backward compatible
2. ✅ Migration can be incremental and gradual
3. ✅ Minimal code changes needed (often just import + parameter)
4. ✅ Test thoroughly before production deployment
5. ✅ Start with simple cases, scale up gradually

**Next steps:**

1. Choose migration strategy (quick win, controlled, or feature flag)
2. Update one code path as proof of concept
3. Test and benchmark
4. Gradually roll out to more code paths
5. Monitor and optimize

Happy parallelizing! 🚀
