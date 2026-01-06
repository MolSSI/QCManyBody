# ParallelConfig API Reference

The `ParallelConfig` class provides configuration options for parallel execution of many-body calculations.

## Class Definition

```python
@dataclass
class ParallelConfig:
    """Configuration for parallel execution parameters.

    This dataclass encapsulates all configuration options for the
    ParallelManyBodyExecutor, including execution mode, resource limits,
    and QCEngine integration settings.
    """
```

## Configuration Parameters

### Execution Control

#### `max_workers: int = 4`
Maximum number of parallel workers for fragment execution.

**Default:** `4`
**Valid Range:** `1` to system CPU count
**Notes:**
- More workers can improve performance but increase memory usage
- Consider memory limits when setting high worker counts
- Start with CPU core count and adjust based on performance

**Example:**
```python
# Conservative setting for limited memory
config = ParallelConfig(max_workers=2)

# Aggressive setting for high-memory systems
config = ParallelConfig(max_workers=8)
```

#### `execution_mode: str = "multiprocessing"`
Parallel execution strategy.

**Default:** `"multiprocessing"`
**Valid Values:**
- `"serial"`: Sequential execution (debugging/validation)
- `"threading"`: Thread-based parallelism (good for I/O-bound tasks)
- `"multiprocessing"`: Process-based parallelism (CPU-intensive tasks)

**Guidelines:**
- **Serial**: Use for debugging or validation against sequential execution
- **Threading**: Often best for QC calculations (I/O-bound with QCEngine)
- **Multiprocessing**: Best for CPU-intensive calculations without external programs

**Example:**
```python
# For debugging/validation
config = ParallelConfig(execution_mode="serial")

# For typical QC calculations
config = ParallelConfig(execution_mode="threading")

# For pure CPU calculations
config = ParallelConfig(execution_mode="multiprocessing")
```

### Resource Management

#### `memory_limit_mb: int = 1000`
Memory limit per worker in megabytes.

**Default:** `1000` MB (1 GB)
**Valid Range:** `100` to system memory / max_workers
**Notes:**
- Passed to QCEngine for memory management
- Should account for QC program memory requirements
- Total memory usage ≈ max_workers × memory_limit_mb

**Example:**
```python
# Conservative for small molecules
config = ParallelConfig(memory_limit_mb=500)

# Generous for large molecules or expensive methods
config = ParallelConfig(memory_limit_mb=4000)
```

#### `timeout_seconds: int = 3600`
Timeout for individual fragment calculations in seconds.

**Default:** `3600` (1 hour)
**Valid Range:** `60` to unlimited
**Notes:**
- Prevents hanging calculations from blocking execution
- Should be set based on expected calculation complexity
- Failed timeouts raise RuntimeError

**Example:**
```python
# Quick calculations (HF/small basis)
config = ParallelConfig(timeout_seconds=300)  # 5 minutes

# Expensive calculations (CCSD(T)/large basis)
config = ParallelConfig(timeout_seconds=7200)  # 2 hours
```

### QCEngine Integration

#### `use_qcengine: bool = True`
Enable or disable QCEngine integration for real quantum chemistry calculations.

**Default:** `True`
**Notes:**
- When `True`: Uses QCEngine to execute real QC calculations
- When `False`: Uses placeholder calculations for testing/validation
- Requires QCEngine and quantum chemistry programs when enabled

**Example:**
```python
# For real calculations
config = ParallelConfig(use_qcengine=True)

# For testing/validation without QC programs
config = ParallelConfig(use_qcengine=False)
```

#### `qc_program: str = "psi4"`
Quantum chemistry program to use via QCEngine.

**Default:** `"psi4"`
**Valid Values:** Any QCEngine-supported program (`"psi4"`, `"nwchem"`, `"cfour"`, etc.)
**Notes:**
- Program must be installed and available to QCEngine
- Different programs may have different performance characteristics
- Program-specific options can be set via `qcengine_config`

**Example:**
```python
# Use Psi4
config = ParallelConfig(qc_program="psi4")

# Use NWChem
config = ParallelConfig(qc_program="nwchem")

# Use CFOUR
config = ParallelConfig(qc_program="cfour")
```

#### `basis_set: str = "sto-3g"`
Default basis set for quantum chemistry calculations.

**Default:** `"sto-3g"`
**Notes:**
- Used when QCEngine integration is disabled (placeholder mode)
- For real calculations, basis set is determined by ManyBodyCore configuration
- Should match the complexity of your intended calculations

#### `qcengine_config: Dict[str, Any] = None`
Additional configuration for QCEngine execution.

**Default:** `None` (empty dictionary)
**Structure:**
```python
{
    "keywords": {},      # QC program-specific keywords
    "protocols": {},     # QCEngine protocols
    "task_config": {}    # Additional task configuration
}
```

**Example:**
```python
config = ParallelConfig(
    qcengine_config={
        "keywords": {
            "scf_type": "df",           # Density fitting for SCF
            "mp2_type": "df",           # Density fitting for MP2
            "freeze_core": True         # Freeze core orbitals
        },
        "protocols": {
            "stdout": False             # Suppress program output
        },
        "task_config": {
            "scratch_directory": "/tmp" # Custom scratch directory
        }
    }
)
```

## Validation and Error Handling

### Constructor Validation

The `__post_init__` method validates configuration parameters:

```python
def __post_init__(self):
    if self.qcengine_config is None:
        self.qcengine_config = {}

    if self.use_qcengine and not HAS_QCENGINE:
        raise ImportError(
            "QCEngine is required for parallel execution but not available. "
            "Install with: pip install qcengine"
        )

    if self.execution_mode not in ["multiprocessing", "threading", "serial"]:
        raise ValueError(f"Invalid execution_mode: {self.execution_mode}")

    if self.max_workers < 1:
        raise ValueError(f"max_workers must be >= 1, got {self.max_workers}")
```

### Common Validation Errors

**ImportError**: Raised when QCEngine is required but not available
```python
# Fix by installing QCEngine
pip install qcengine
```

**ValueError**: Raised for invalid parameter values
```python
# Invalid execution mode
config = ParallelConfig(execution_mode="invalid")  # ❌ ValueError

# Invalid worker count
config = ParallelConfig(max_workers=0)  # ❌ ValueError
```

## Configuration Patterns

### Development and Testing

```python
# For development/testing without QC programs
dev_config = ParallelConfig(
    max_workers=1,
    execution_mode="serial",
    use_qcengine=False,
    timeout_seconds=60
)
```

### Production Calculations

```python
# For production calculations with Psi4
prod_config = ParallelConfig(
    max_workers=4,
    execution_mode="threading",
    use_qcengine=True,
    qc_program="psi4",
    memory_limit_mb=2000,
    timeout_seconds=3600,
    qcengine_config={
        "keywords": {"scf_type": "df"},
        "protocols": {"stdout": False}
    }
)
```

### HPC Environment

```python
# For HPC clusters with many cores
hpc_config = ParallelConfig(
    max_workers=16,
    execution_mode="multiprocessing",
    use_qcengine=True,
    memory_limit_mb=4000,
    timeout_seconds=7200,
    qcengine_config={
        "task_config": {
            "scratch_directory": "/scratch/$SLURM_JOB_ID"
        }
    }
)
```

### Memory-Constrained Systems

```python
# For systems with limited memory
memory_config = ParallelConfig(
    max_workers=2,
    execution_mode="threading",
    memory_limit_mb=500,
    timeout_seconds=1800
)
```

## Best Practices

### Worker Count Selection

1. **Start Conservative**: Begin with 2-4 workers
2. **Monitor Memory**: Ensure total memory usage < system memory
3. **Consider QC Program**: Some programs scale better than others
4. **Test and Adjust**: Profile with your specific calculations

### Memory Management

1. **Per-Worker Limits**: Set `memory_limit_mb` based on largest expected fragment
2. **Total Memory Budget**: Ensure `max_workers × memory_limit_mb < system_memory`
3. **QC Program Overhead**: Account for program-specific memory requirements
4. **Monitor Usage**: Check actual memory consumption during execution

### Execution Mode Selection

1. **Threading for QCEngine**: Usually best for QC calculations
2. **Multiprocessing for CPU**: Best for pure computational tasks
3. **Serial for Debugging**: Always test with serial mode first
4. **Consider I/O**: Threading often better for I/O-bound tasks

### Timeout Configuration

1. **Method Complexity**: Longer timeouts for expensive methods (CCSD(T), etc.)
2. **System Size**: Larger molecules need longer timeouts
3. **Basis Set Size**: Larger basis sets increase calculation time
4. **Safety Margin**: Set 2-3x expected calculation time

---

For usage examples and patterns, see the [Usage Guide](../usage/README.md).