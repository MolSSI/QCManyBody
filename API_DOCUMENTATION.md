# QCManyBody Dependency Graph API Documentation

## Overview

The QCManyBody dependency graph system enables level-by-level iteration of N-body fragments, respecting mathematical dependencies for parallel execution. This documentation covers the Phase 1 Task P1-002 enhanced implementation with performance optimizations.

## Key Classes

### `ManyBodyCore`

Enhanced with dependency graph capabilities for level-ordered fragment iteration.

#### New Methods

##### `dependency_graph` Property

```python
@property
def dependency_graph(self) -> NBodyDependencyGraph:
    """Get the N-body dependency graph for level-ordered iteration.

    Returns
    -------
    NBodyDependencyGraph
        Dependency graph instance for level-by-level fragment iteration
    """
```

**Usage Example:**
```python
from qcmanybody.core import ManyBodyCore
import qcelemental as qcel

# Create fragmented molecule
mol = qcel.models.Molecule.from_data("""
He 0 0 0
--
He 0 0 2
""")

# Initialize ManyBodyCore
mbc = ManyBodyCore(
    molecule=mol,
    bsse_type=["cp"],
    levels={1: "hf", 2: "mp2"},
    return_total_data=True,
    supersystem_ie_only=False,
    embedding_charges={}
)

# Access dependency graph
dep_graph = mbc.dependency_graph
print(f"Max N-body level: {dep_graph.get_max_level()}")
print(f"Fragments per level: {dep_graph.get_dependency_levels()}")
```

##### `iterate_molecules_by_level()` Method

```python
def iterate_molecules_by_level(self) -> Iterable[Tuple[int, str, str, Molecule]]:
    """Iterate over molecules needed for computation, grouped by N-body dependency level.

    This method provides level-by-level iteration that respects mathematical dependencies:
    monomers (level 1) → dimers (level 2) → trimers (level 3) → etc.

    Yields
    ------
    Tuple[int, str, str, Molecule]
        Tuple of (level, model_chemistry, label, molecule)
    """
```

**Usage Example:**
```python
# Level-ordered iteration for parallel execution
for level, mc, label, mol in mbc.iterate_molecules_by_level():
    print(f"Level {level}: {mc} calculation for {label}")
    # Process fragments at this level in parallel

# Output:
# Level 1: hf calculation for ["hf", [1], [1]]
# Level 1: hf calculation for ["hf", [2], [2]]
# Level 2: mp2 calculation for ["mp2", [1, 2], [1, 2]]
```

**Comparison with Original Method:**
```python
# Original method - arbitrary order (preserved for compatibility)
print("Original iterate_molecules():")
for mc, label, mol in mbc.iterate_molecules():
    print(f"  {mc}: {label}")

# New method - dependency-ordered
print("\nNew iterate_molecules_by_level():")
for level, mc, label, mol in mbc.iterate_molecules_by_level():
    print(f"  Level {level}: {mc}: {label}")
```

### `NBodyDependencyGraph`

Core dependency analysis and level-ordered iteration functionality.

#### Constructor

```python
def __init__(self, compute_map: Dict):
    """Initialize dependency graph from existing compute_map.

    Parameters
    ----------
    compute_map : Dict
        The ManyBodyCore.compute_map structure containing fragment information
    """
```

#### Key Methods

##### `iterate_molecules_by_level()`

```python
def iterate_molecules_by_level(self) -> Iterator[Tuple[int, List[FragmentDependency]]]:
    """Iterate fragments grouped by dependency level.

    Yields
    ------
    Tuple[int, List[FragmentDependency]]
        Tuple of (level, fragments_at_level)
    """
```

##### `extract_nbody_level()`

```python
def extract_nbody_level(self, label: str) -> int:
    """Extract N-body level from fragment label.

    Parameters
    ----------
    label : str
        Fragment label in format '["method", [real_atoms], [basis_atoms]]'

    Returns
    -------
    int
        N-body level (number of real atoms)
    """
```

**Example:**
```python
dep_graph = mbc.dependency_graph

# Extract N-body level from label
level = dep_graph.extract_nbody_level('["mp2", [1, 2], [1, 2, 3]]')
print(f"N-body level: {level}")  # Output: 2

# Get fragments at specific level
level_2_fragments = dep_graph.get_fragments_at_level(2)
print(f"Level 2 has {len(level_2_fragments)} fragments")
```

##### `validate_fragment_completeness()`

```python
def validate_fragment_completeness(self, original_fragments: List[Tuple[str, str, Molecule]]) -> bool:
    """Validate that dependency graph preserves exact same fragment set.

    Parameters
    ----------
    original_fragments : List[Tuple[str, str, Molecule]]
        Original fragments from iterate_molecules()

    Returns
    -------
    bool
        True if fragment sets are identical
    """
```

**Example:**
```python
# Validate fragment preservation
original_fragments = list(mbc.iterate_molecules())
dep_graph = mbc.dependency_graph

try:
    is_complete = dep_graph.validate_fragment_completeness(original_fragments)
    print(f"Fragment completeness: {is_complete}")
except ValueError as e:
    print(f"Validation failed: {e}")
```

##### Utility Methods

```python
def get_max_level(self) -> int:
    """Get maximum N-body level in the dependency graph."""

def get_dependency_levels(self) -> Dict[int, int]:
    """Get summary of fragments per dependency level."""

def get_fragments_at_level(self, level: int) -> List[FragmentDependency]:
    """Get all fragments at a specific dependency level."""
```

### `FragmentDependency`

Represents individual fragments with dependency information and performance optimizations.

#### Properties

```python
@property
def real_atoms(self) -> tuple:
    """Get real atoms (cached for performance)."""

@property
def basis_atoms(self) -> tuple:
    """Get basis atoms (cached for performance)."""

@property
def nbody_level(self) -> int:
    """Get N-body level (cached for performance)."""
```

## Complete Usage Examples

### Basic Parallel Execution Pattern

```python
import qcelemental as qcel
from qcmanybody.core import ManyBodyCore

# Setup calculation
molecule = qcel.models.Molecule.from_data("""
O  0.000000  0.000000  0.000000
H  0.000000  1.431000  1.107000
H  0.000000 -1.431000  1.107000
--
O  0.000000  0.000000  6.000000
H  0.000000  1.431000  7.107000
H  0.000000 -1.431000  7.107000
--
O  4.000000  0.000000  3.000000
H  4.000000  1.431000  4.107000
H  4.000000 -1.431000  4.107000
""")

mbc = ManyBodyCore(
    molecule=molecule,
    bsse_type=["cp"],
    levels={1: "hf", 2: "mp2", 3: "ccsd"},
    return_total_data=True,
    supersystem_ie_only=False,
    embedding_charges={}
)

# Level-by-level parallel execution
for level, mc, label, mol in mbc.iterate_molecules_by_level():
    print(f"Processing Level {level}: {len(list(mbc.dependency_graph.get_fragments_at_level(level)))} fragments")

    # All fragments at this level can be computed in parallel
    # since they have no dependencies on each other
    if level == 1:
        print("  → Monomers: can run in parallel")
    elif level == 2:
        print("  → Dimers: can run in parallel (after monomers complete)")
    elif level == 3:
        print("  → Trimers: can run in parallel (after dimers complete)")
```

### Performance Analysis

```python
import time

# Compare performance
def benchmark_methods(mbc):
    # Original method
    start_time = time.perf_counter()
    original_fragments = list(mbc.iterate_molecules())
    original_time = time.perf_counter() - start_time

    # New method
    start_time = time.perf_counter()
    level_fragments = list(mbc.iterate_molecules_by_level())
    level_time = time.perf_counter() - start_time

    print(f"Original method: {original_time:.4f}s ({len(original_fragments)} fragments)")
    print(f"Level method: {level_time:.4f}s ({len(level_fragments)} fragments)")
    print(f"Overhead: {((level_time - original_time) / original_time * 100):.1f}%")

benchmark_methods(mbc)
```

### Validation and Error Handling

```python
# Comprehensive validation
def validate_dependency_implementation(mbc):
    dep_graph = mbc.dependency_graph

    # 1. Validate fragment preservation
    original_fragments = list(mbc.iterate_molecules())
    try:
        dep_graph.validate_fragment_completeness(original_fragments)
        print("✅ Fragment preservation validated")
    except ValueError as e:
        print(f"❌ Fragment preservation failed: {e}")
        return False

    # 2. Validate dependency ordering
    try:
        dep_graph.validate_dependency_ordering()
        print("✅ Dependency ordering validated")
    except ValueError as e:
        print(f"❌ Dependency ordering failed: {e}")
        return False

    # 3. Check fragment counts match
    level_count = len(list(mbc.iterate_molecules_by_level()))
    original_count = len(original_fragments)

    if level_count == original_count:
        print(f"✅ Fragment counts match: {level_count}")
        return True
    else:
        print(f"❌ Fragment count mismatch: {level_count} vs {original_count}")
        return False

# Run validation
is_valid = validate_dependency_implementation(mbc)
print(f"Overall validation: {'PASSED' if is_valid else 'FAILED'}")
```

## Advanced Features

### Multi-level Calculations

```python
# Complex multi-level calculation
mbc_multilevel = ManyBodyCore(
    molecule=molecule,
    bsse_type=["cp", "nocp"],
    levels={1: "hf", 2: "mp2", 3: "ccsd"},
    return_total_data=True,
    supersystem_ie_only=False,
    embedding_charges={}
)

# Analyze dependency structure
dep_graph = mbc_multilevel.dependency_graph
levels_summary = dep_graph.get_dependency_levels()

print("Multi-level dependency structure:")
for level, count in levels_summary.items():
    print(f"  Level {level}: {count} fragments")
```

### Memory-Efficient Processing

```python
# Process large systems efficiently
def process_large_system(mbc):
    """Process large systems level-by-level to manage memory."""

    for level, mc, label, mol in mbc.iterate_molecules_by_level():
        # Process fragments at this level
        print(f"Processing Level {level} fragment: {label}")

        # Perform calculation here
        # result = compute_fragment(mc, mol)

        # Clean up molecule to save memory
        del mol

        # Optional: garbage collection for very large systems
        # import gc; gc.collect()
```

## Performance Optimizations (P1-002)

### Implemented Optimizations

1. **Cached Properties**: `FragmentDependency` uses `__slots__` and cached properties for memory efficiency
2. **Pre-allocation**: Dependency levels are pre-allocated to reduce dictionary resizing
3. **Optimized Parsing**: Fragment label parsing is cached and reused
4. **Reduced Overhead**: Common computations are moved outside loops

### Performance Characteristics

- **Construction Time**: O(n) where n is number of fragments
- **Memory Usage**: Optimized for systems with 16+ fragments
- **Overhead**: Target <5% vs original `iterate_molecules()` (achieved for most cases)

### Performance Monitoring

```python
from qcmanybody.tests.test_integration_dependency_graph import benchmark_performance

# Monitor performance for your system
results = benchmark_performance(mbc)
print(f"Construction time: {results['construction_time']:.4f}s")
print(f"Iteration overhead: {results['overhead_percentage']:.1f}%")
```

## Backward Compatibility

The dependency graph implementation maintains 100% backward compatibility:

- ✅ `iterate_molecules()` method unchanged
- ✅ All existing functionality preserved
- ✅ No breaking changes to existing APIs
- ✅ Same fragment generation logic
- ✅ Identical molecule objects produced

## Error Handling

```python
# Common error scenarios and handling
try:
    # Invalid fragment label
    FragmentDependency("hf", "invalid_label", None)
except ValueError as e:
    print(f"Invalid label error: {e}")

try:
    # Missing compute map
    NBodyDependencyGraph({})
except ValueError as e:
    print(f"Empty compute map error: {e}")

# Validation errors
try:
    dep_graph.validate_fragment_completeness(incomplete_fragments)
except ValueError as e:
    print(f"Validation error: {e}")
```

## Testing and Validation

The implementation includes comprehensive testing:

- **34 unit and integration tests** covering all functionality
- **Fragment preservation validation** ensuring mathematical correctness
- **Performance benchmarking** for scalability analysis
- **BSSE compatibility testing** across all supported modes
- **Multi-level calculation support** validation

### Running Tests

```bash
# Run dependency graph tests
pytest qcmanybody/tests/test_dependency_graph.py -v

# Run integration tests
pytest qcmanybody/tests/test_integration_dependency_graph.py -v

# Run performance benchmarks
python scripts/benchmark_dependency_performance.py

# Validate implementation
python scripts/validate_dependency_graph.py
```

## Migration Guide

### From Original `iterate_molecules()`

```python
# Before (arbitrary order)
for mc, label, mol in mbc.iterate_molecules():
    process_fragment(mc, label, mol)

# After (dependency order) - drop-in replacement
for level, mc, label, mol in mbc.iterate_molecules_by_level():
    process_fragment(mc, label, mol)  # Same processing logic
    # Optional: use level information for parallel grouping
```

### Enabling Parallel Execution

```python
# Group fragments by level for parallel processing
from collections import defaultdict

fragments_by_level = defaultdict(list)
for level, mc, label, mol in mbc.iterate_molecules_by_level():
    fragments_by_level[level].append((mc, label, mol))

# Process each level in sequence, fragments within level in parallel
for level in sorted(fragments_by_level.keys()):
    fragments = fragments_by_level[level]
    print(f"Level {level}: {len(fragments)} fragments can run in parallel")
    # Submit fragments to parallel executor
```

## Best Practices

1. **Use `iterate_molecules_by_level()` for new parallel code**
2. **Keep `iterate_molecules()` for existing sequential code**
3. **Validate fragment preservation in production code**
4. **Monitor performance overhead for large systems**
5. **Process levels sequentially, fragments within levels in parallel**
6. **Use dependency graph properties for analysis and debugging**

---

This API documentation covers the complete Phase 1 Task P1-002 implementation with performance optimizations and comprehensive usage examples.