# Quick Start Guide

This guide gets you up and running with parallel execution in QCManyBody in just a few minutes.

## 📦 Installation Requirements

### Basic Requirements
```bash
# Core QCManyBody (already installed)
pip install qcelemental numpy pydantic

# For parallel execution (built-in)
# No additional packages needed for basic functionality
```

### For Real Quantum Chemistry
```bash
# Install QCEngine
pip install qcengine

# Install a quantum chemistry program (choose one)
conda install -c conda-forge psi4        # Psi4
conda install -c conda-forge nwchem      # NWChem
# or other QCEngine-supported programs
```

## 🚀 Your First Parallel Calculation

### Step 1: Import Required Modules

```python
from qcmanybody import ManyBodyCore, BsseEnum
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
from qcelemental.models import Molecule
```

### Step 2: Create a Molecular System

```python
# Simple water dimer
water_dimer = Molecule.from_data("""
O  0.0000  0.0000  0.0000
H  0.7570  0.5860  0.0000
H -0.7570  0.5860  0.0000
--
O  3.0000  0.0000  0.0000
H  3.7570  0.5860  0.0000
H  2.2430  0.5860  0.0000
""")

print(f"Created system with {len(water_dimer.fragments)} fragments")
```

### Step 3: Set Up Many-Body Calculation

```python
# Configure the many-body calculation
core = ManyBodyCore(
    molecule=water_dimer,
    bsse_type=[BsseEnum.nocp],  # No counterpoise correction
    levels={1: "hf", 2: "hf"},  # HF for both monomers and dimers
    return_total_data=False,
    supersystem_ie_only=False,
    embedding_charges={}
)

print(f"Many-body levels: {core.levels}")
```

### Step 4: Configure Parallel Execution

```python
# Start with a simple configuration
config = ParallelConfig(
    max_workers=2,              # Use 2 parallel workers
    execution_mode="threading", # Threading mode (good for QC)
    use_qcengine=False,        # Start with placeholder calculations
    basis_set="sto-3g"         # Simple basis set
)

print(f"Parallel config: {config.max_workers} workers, {config.execution_mode} mode")
```

### Step 5: Execute the Calculation

```python
# Create executor and run calculation
executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()

print(f"Calculation completed: {len(results)} fragments calculated")
```

### Step 6: Examine Results

```python
# Print results
print("\nFragment Results:")
for label, result in results.items():
    print(f"  {label}: {result.return_result:.8f} hartree")

# Get execution statistics
stats = executor.get_execution_statistics()
print(f"\nExecution Statistics:")
print(f"  Total fragments: {stats['total_fragments']}")
print(f"  Levels executed: {stats['levels_executed']}")
print(f"  Execution time: {stats['parallel_time']:.3f}s")
print(f"  Estimated speedup: {stats['speedup_factor']:.2f}x")
```

## ✅ Complete First Example

Here's the complete code for your first parallel calculation:

```python
from qcmanybody import ManyBodyCore, BsseEnum
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
from qcelemental.models import Molecule

# Create water dimer
water_dimer = Molecule.from_data("""
O  0.0000  0.0000  0.0000
H  0.7570  0.5860  0.0000
H -0.7570  0.5860  0.0000
--
O  3.0000  0.0000  0.0000
H  3.7570  0.5860  0.0000
H  2.2430  0.5860  0.0000
""")

# Set up many-body calculation
core = ManyBodyCore(
    molecule=water_dimer,
    bsse_type=[BsseEnum.nocp],
    levels={1: "hf", 2: "hf"},
    return_total_data=False,
    supersystem_ie_only=False,
    embedding_charges={}
)

# Configure parallel execution
config = ParallelConfig(
    max_workers=2,
    execution_mode="threading",
    use_qcengine=False  # Placeholder mode for quick testing
)

# Execute calculation
executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()

# Print results
print("🚀 Parallel Many-Body Calculation Results")
print("=" * 45)
for label, result in results.items():
    print(f"{label}: {result.return_result:.8f} hartree")

stats = executor.get_execution_statistics()
print(f"\n📊 Statistics:")
print(f"Fragments: {stats['total_fragments']}, Time: {stats['parallel_time']:.3f}s")
```

**Expected Output:**
```
🚀 Parallel Many-Body Calculation Results
=============================================
["hf", [1], [1]]: -3.000000 hartree
["hf", [2], [2]]: -3.000000 hartree
["hf", [1, 2], [1, 2]]: -6.000000 hartree

📊 Statistics:
Fragments: 3, Time: 0.023s
```

## 🧪 Testing with Real Quantum Chemistry

Once you have QCEngine and a quantum chemistry program installed:

```python
# Enable real QC calculations
config = ParallelConfig(
    max_workers=2,
    execution_mode="threading",
    use_qcengine=True,      # Enable QCEngine
    qc_program="psi4",      # Use Psi4
    basis_set="sto-3g",
    memory_limit_mb=1000,
    timeout_seconds=300
)

# Execute with real quantum chemistry
executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()

print("Real QC Results:")
for label, result in results.items():
    print(f"{label}: {result.return_result:.8f} hartree")
```

## ⚙️ Multiprocessing start methods

QCManyBody now supports Python's ``spawn`` start method for multiprocessing,
which is the default on Windows and macOS and available on modern Linux
distributions. The parallel executor ships with a lightweight initializer that
hydrates worker processes with the serialized ``ParallelConfig`` so fragment
payloads remain picklable under both ``fork`` and ``spawn``.

- **Testing tip:** use ``ParallelConfig(use_qcengine=False)`` to exercise the
    placeholder workflow without requiring Psi4. The regression test
    ``qcmanybody/tests/test_multiprocessing_serialization.py`` loads the
    real ``test_water4_mbe4_multiprocessing.py`` example and validates spawn-mode
    execution against a serial baseline.
- **Production tip:** set the start method explicitly early in your program via
    ``multiprocessing.set_start_method("spawn", force=True)`` when targeting
    heterogeneous platforms, and ensure any heavy dependencies are imported
    inside worker initializers rather than at module scope.

When ``use_qcengine=True``, the executor continues to require the relevant
quantum-chemistry backends (Psi4, NWChem, etc.) to be installed in the runtime
environment.

## 🔍 Validation Example

Verify that parallel execution gives identical results to sequential:

```python
# Sequential execution
config_seq = ParallelConfig(
    execution_mode="serial",
    max_workers=1,
    use_qcengine=False
)

executor_seq = ParallelManyBodyExecutor(core, config_seq)
results_seq = executor_seq.execute_full_calculation()

# Parallel execution
config_par = ParallelConfig(
    execution_mode="threading",
    max_workers=2,
    use_qcengine=False
)

executor_par = ParallelManyBodyExecutor(core, config_par)
results_par = executor_par.execute_full_calculation()

# Validate correctness
is_correct = executor_par.validate_parallel_correctness(
    results_par, results_seq, tolerance=1e-12
)

print(f"✅ Parallel correctness: {'PASSED' if is_correct else 'FAILED'}")
```

## 🎯 Next Steps

Now that you have parallel execution working:

1. **Try Different Systems**: Test with water trimers, tetramers, etc.
2. **Experiment with Settings**: Adjust `max_workers`, `execution_mode`
3. **Enable Real QC**: Install Psi4 and try `use_qcengine=True`
4. **Monitor Performance**: Use `get_execution_statistics()` to track speedup
5. **Explore Advanced Features**: Check out the [Advanced Usage](advanced-usage.md) guide

## 🔧 Troubleshooting Quick Start

### Common Issues

**ImportError: No module named 'qcmanybody.parallel'**
```python
# Make sure you're using the latest version with parallel support
# Check that qcmanybody/parallel.py exists in your installation
```

**RuntimeError: P1-002 dependency graph foundation required**
```python
# This should not happen with current version
# If it does, check that ManyBodyCore has iterate_molecules_by_level method
```

**Memory or timeout errors with QCEngine**
```python
# Start with placeholder mode (use_qcengine=False)
# Then gradually enable QCEngine with small systems
```

### Getting Help

- **Basic Usage**: Continue to [Basic Examples](basic-examples.md)
- **Advanced Features**: See [Advanced Usage](advanced-usage.md)
- **Production Setup**: Check [Production Workflows](production-workflows.md)
- **Issues**: Review [Troubleshooting](troubleshooting.md)

---

🎉 **Congratulations!** You've successfully run your first parallel many-body calculation with QCManyBody!