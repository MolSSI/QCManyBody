# Usage Guide

This guide provides practical examples and usage patterns for the QCManyBody Parallel Execution system.

## 📚 Usage Documentation

### Getting Started
- **[Quick Start](quick-start.md)** - Basic usage patterns
- **[Basic Examples](basic-examples.md)** - Simple examples for common use cases
- **[Advanced Usage](advanced-usage.md)** - Complex scenarios and optimization

### Practical Applications
- **[Production Workflows](production-workflows.md)** - Real-world calculation setups
- **[HPC Deployment](hpc-deployment.md)** - High-performance computing environments
- **[Troubleshooting](troubleshooting.md)** - Common issues and solutions

## 🚀 Quick Start

### Basic Parallel Execution

```python
from qcmanybody import ManyBodyCore, BsseEnum
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
from qcelemental.models import Molecule

# 1. Create your molecular system
water_dimer = Molecule.from_data("""
O  0.0000  0.0000  0.0000
H  0.7570  0.5860  0.0000
H -0.7570  0.5860  0.0000
--
O  3.0000  0.0000  0.0000
H  3.7570  0.5860  0.0000
H  2.2430  0.5860  0.0000
""")

# 2. Set up many-body calculation
core = ManyBodyCore(
    molecule=water_dimer,
    bsse_type=[BsseEnum.nocp],
    levels={1: "hf", 2: "mp2"},
    return_total_data=False,
    supersystem_ie_only=False,
    embedding_charges={}
)

# 3. Configure parallel execution
config = ParallelConfig(
    max_workers=4,
    execution_mode="threading",
    use_qcengine=True,
    qc_program="psi4"
)

# 4. Execute with parallelization
executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()

# 5. Access results
for label, result in results.items():
    print(f"{label}: {result.return_result:.8f} hartree")
```

### Key Benefits

- **🚀 Performance**: Significant speedup for multi-fragment systems
- **🔬 Precision**: Maintains 1e-12 quantum chemistry accuracy
- **🔗 Dependencies**: Respects mathematical N-body dependencies
- **⚙️ Flexibility**: Multiple execution modes and configurations
- **🧪 Validation**: Built-in correctness verification

## 📋 Common Usage Patterns

### 1. Development and Testing

```python
# Use placeholder execution for fast testing
config = ParallelConfig(
    max_workers=1,
    execution_mode="serial",
    use_qcengine=False  # Placeholder calculations
)

executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()

# Validate against known reference
is_correct = executor.validate_parallel_correctness(
    results, reference_results, tolerance=1e-12
)
```

### 2. Production Calculations

```python
# Optimized for real quantum chemistry
config = ParallelConfig(
    max_workers=8,
    execution_mode="threading",
    use_qcengine=True,
    qc_program="psi4",
    memory_limit_mb=2000,
    timeout_seconds=3600,
    qcengine_config={
        "keywords": {"scf_type": "df", "mp2_type": "df"},
        "protocols": {"stdout": False}
    }
)

executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()

# Monitor performance
stats = executor.get_execution_statistics()
print(f"Speedup: {stats['speedup_factor']:.2f}x")
```

### 3. Validation and Verification

```python
# Compare parallel vs sequential execution
configs = {
    "sequential": ParallelConfig(execution_mode="serial", max_workers=1),
    "parallel": ParallelConfig(execution_mode="threading", max_workers=4)
}

results = {}
for name, config in configs.items():
    executor = ParallelManyBodyExecutor(core, config)
    results[name] = executor.execute_full_calculation()

# Validate mathematical correctness
executor.validate_parallel_correctness(
    results["parallel"], results["sequential"], tolerance=1e-12
)
```

## 🎯 Use Case Examples

### Small Molecules (2-3 fragments)
```python
# Conservative settings for small systems
config = ParallelConfig(
    max_workers=2,
    execution_mode="threading",
    memory_limit_mb=500,
    timeout_seconds=300
)
```

### Medium Molecules (4-6 fragments)
```python
# Balanced settings for medium systems
config = ParallelConfig(
    max_workers=4,
    execution_mode="threading",
    memory_limit_mb=1000,
    timeout_seconds=1800
)
```

### Large Molecules (7+ fragments)
```python
# Aggressive settings for large systems
config = ParallelConfig(
    max_workers=8,
    execution_mode="threading",
    memory_limit_mb=2000,
    timeout_seconds=3600
)
```

## ⚡ Performance Tips

### 1. Worker Count Optimization
- Start with CPU core count
- Monitor memory usage
- Adjust based on actual performance

### 2. Memory Management
- Set `memory_limit_mb` based on largest fragment
- Ensure total memory < system memory
- Consider QC program memory overhead

### 3. Execution Mode Selection
- **Threading**: Best for QCEngine integration
- **Multiprocessing**: Best for CPU-intensive pure calculations
- **Serial**: For debugging and validation

### 4. QCEngine Configuration
- Use density fitting (`scf_type="df"`) for speed
- Suppress output (`stdout=False`) for clean execution
- Set appropriate scratch directories

## 🔧 Error Handling

### Common Issues and Solutions

**Memory Errors:**
```python
# Reduce memory usage
config = ParallelConfig(
    max_workers=2,  # Fewer workers
    memory_limit_mb=500  # Less memory per worker
)
```

**Timeout Errors:**
```python
# Increase timeout
config = ParallelConfig(
    timeout_seconds=7200  # 2 hours
)
```

**QCEngine Errors:**
```python
# Use placeholder mode for testing
config = ParallelConfig(
    use_qcengine=False  # Skip QCEngine
)
```

## 📊 Monitoring and Statistics

### Execution Statistics

```python
# Get detailed performance metrics
stats = executor.get_execution_statistics()

print(f"Total fragments: {stats['total_fragments']}")
print(f"Levels executed: {stats['levels_executed']}")
print(f"Parallel time: {stats['parallel_time']:.3f}s")
print(f"Estimated speedup: {stats['speedup_factor']:.2f}x")
```

### Progress Monitoring

```python
import logging

# Enable debug logging for detailed progress
logging.basicConfig(level=logging.DEBUG)

# Execute with detailed logging
executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()
```

## 🔗 Integration Examples

### With Existing Workflows

```python
def run_parallel_manybody(molecule, method_levels, **kwargs):
    """Convenience function for parallel many-body calculations."""

    # Set up ManyBodyCore
    core = ManyBodyCore(
        molecule=molecule,
        bsse_type=[BsseEnum.nocp],
        levels=method_levels,
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    # Configure parallel execution
    config = ParallelConfig(**kwargs)

    # Execute and return results
    executor = ParallelManyBodyExecutor(core, config)
    return executor.execute_full_calculation()

# Use in workflows
results = run_parallel_manybody(
    molecule=water_trimer,
    method_levels={1: "hf", 2: "mp2", 3: "ccsd(t)"},
    max_workers=8,
    execution_mode="threading"
)
```

---

For more detailed examples, see the individual usage documentation pages.