# API Reference

This section provides comprehensive documentation for the QCManyBody Parallel Execution API.

## 📚 API Documentation Structure

### Core Classes
- **[ParallelManyBodyExecutor](parallel-executor.md)** - Main parallel execution engine
- **[ParallelConfig](parallel-config.md)** - Configuration for parallel execution

### Validation & Testing
- **[Validation Framework](validation-framework.md)** - Ultra-strict correctness validation
- **[Testing Utilities](testing-utilities.md)** - Testing and benchmarking tools

### Examples & Integration
- **[Code Examples](examples.md)** - Practical usage examples
- **[QCEngine Integration](qcengine-integration.md)** - Quantum chemistry program integration

## 🚀 Quick Reference

### Basic Usage

```python
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig

# Configure parallel execution
config = ParallelConfig(
    max_workers=4,
    execution_mode="threading",
    use_qcengine=True,
    qc_program="psi4"
)

# Create and execute
executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()
```

### Key Methods

| Method | Description | Returns |
|--------|-------------|---------|
| `execute_full_calculation()` | Execute complete many-body calculation | `Dict[str, AtomicResult]` |
| `execute_level_parallel()` | Execute fragments at specific level | `Dict[str, AtomicResult]` |
| `get_execution_statistics()` | Get performance metrics | `Dict[str, Union[int, float]]` |
| `validate_parallel_correctness()` | Validate against reference | `bool` |

### Configuration Options

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `max_workers` | `int` | Maximum parallel workers | `4` |
| `execution_mode` | `str` | Execution mode: serial/threading/multiprocessing | `"multiprocessing"` |
| `use_qcengine` | `bool` | Enable QCEngine integration | `True` |
| `qc_program` | `str` | Quantum chemistry program | `"psi4"` |
| `memory_limit_mb` | `int` | Memory limit per worker (MB) | `1000` |
| `timeout_seconds` | `int` | Timeout for calculations (seconds) | `3600` |

## 📊 Return Types

### AtomicResult
Standard QCElemental result format:
```python
{
    "driver": "energy",
    "model": {"method": "hf", "basis": "sto-3g"},
    "molecule": Molecule(...),
    "return_result": -1.1167383286909045,  # Energy in hartree
    "success": True,
    "properties": {},
    "provenance": {"creator": "qcmanybody-parallel", "version": "dev"}
}
```

### Execution Statistics
Performance metrics returned by `get_execution_statistics()`:
```python
{
    "total_fragments": 7,
    "levels_executed": 3,
    "parallel_time": 0.045,
    "sequential_time_estimate": 0.105,
    "speedup_factor": 2.33
}
```

## ⚠️ Important Notes

### Mathematical Correctness
- All parallel results must reproduce sequential results within 1e-12 tolerance
- Fragment execution order within levels does not affect results
- Level execution order is strictly enforced (1 → 2 → 3 → N)

### QCEngine Integration
- Requires QCEngine and quantum chemistry programs (Psi4, NWChem, etc.)
- Memory and core settings are managed automatically
- Program-specific configurations can be passed via `qcengine_config`

### Error Handling
- Fragment failures propagate to level failure
- Level failures stop the entire calculation
- Comprehensive error messages with fragment identification
- Timeout handling for long-running calculations

---

See individual API documentation pages for detailed method signatures, parameters, and examples.