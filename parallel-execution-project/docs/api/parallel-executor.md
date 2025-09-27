# ParallelManyBodyExecutor API Reference

The `ParallelManyBodyExecutor` is the main class for parallel execution of many-body calculations in QCManyBody.

## Class Definition

```python
class ParallelManyBodyExecutor:
    """Parallel execution engine for many-body calculations.

    This class implements level-by-level parallel execution that respects N-body
    dependencies while enabling parallelization within each level. It builds on
    the P1-002 dependency graph foundation to ensure mathematical correctness.
    """
```

## Constructor

### `__init__(core: ManyBodyCore, config: ParallelConfig)`

Initialize the parallel executor.

**Parameters:**
- `core` (`ManyBodyCore`): ManyBodyCore instance with P1-002 dependency graph foundation
- `config` (`ParallelConfig`): Parallel execution configuration

**Raises:**
- `RuntimeError`: If ManyBodyCore missing `iterate_molecules_by_level()` method

**Example:**
```python
from qcmanybody import ManyBodyCore, BsseEnum
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig

core = ManyBodyCore(
    molecule=water_dimer,
    bsse_type=[BsseEnum.nocp],
    levels={1: "hf", 2: "mp2"},
    return_total_data=False,
    supersystem_ie_only=False,
    embedding_charges={}
)

config = ParallelConfig(max_workers=4, execution_mode="threading")
executor = ParallelManyBodyExecutor(core, config)
```

## Main Execution Methods

### `execute_full_calculation() -> Dict[str, AtomicResult]`

Execute complete many-body calculation with level-by-level parallelism.

**Returns:**
- `Dict[str, AtomicResult]`: Dictionary mapping fragment labels to calculation results

**Raises:**
- `RuntimeError`: If parallel execution fails or produces invalid results

**Example:**
```python
results = executor.execute_full_calculation()

# Access results by fragment label
monomer_result = results['["hf", [1], [1]]']
dimer_result = results['["mp2", [1, 2], [1, 2]]']

print(f"Monomer energy: {monomer_result.return_result}")
print(f"Dimer energy: {dimer_result.return_result}")
```

### `execute_level_parallel(level: int, fragments_at_level: List[Tuple]) -> Dict[str, AtomicResult]`

Execute all fragments at a given level in parallel.

**Parameters:**
- `level` (`int`): The N-body dependency level (1 for monomers, 2 for dimers, etc.)
- `fragments_at_level` (`List[Tuple[int, str, str, Molecule]]`): Fragment specifications

**Returns:**
- `Dict[str, AtomicResult]`: Dictionary mapping fragment labels to results

**Raises:**
- `RuntimeError`: If any fragment calculation fails

**Example:**
```python
# Get fragments at level 1 (monomers)
level_1_fragments = []
for level, mc, label, mol in core.iterate_molecules_by_level():
    if level == 1:
        level_1_fragments.append((level, mc, label, mol))

# Execute level 1 in parallel
level_1_results = executor.execute_level_parallel(1, level_1_fragments)
```

## Fragment Execution

### `execute_fragment(fragment_spec: Tuple[int, str, str, Molecule]) -> Tuple[str, AtomicResult]`

Execute a single fragment calculation.

**Parameters:**
- `fragment_spec` (`Tuple[int, str, str, Molecule]`): Fragment specification (level, model_chemistry, label, molecule)

**Returns:**
- `Tuple[str, AtomicResult]`: Tuple of (label, result) for the completed calculation

**Raises:**
- `RuntimeError`: If fragment calculation fails or times out

**Example:**
```python
# Execute single fragment
fragment_spec = (1, "hf", '["hf", [1], [1]]', molecule)
label, result = executor.execute_fragment(fragment_spec)

print(f"Fragment {label}: {result.return_result} hartree")
```

## Statistics and Monitoring

### `get_execution_statistics() -> Dict[str, Union[int, float]]`

Get detailed execution statistics.

**Returns:**
Dictionary containing execution performance metrics:
- `total_fragments` (`int`): Total number of fragments executed
- `levels_executed` (`int`): Number of dependency levels executed
- `parallel_time` (`float`): Total parallel execution time (seconds)
- `sequential_time_estimate` (`float`): Estimated sequential time (seconds)
- `speedup_factor` (`float`): Estimated speedup factor

**Example:**
```python
stats = executor.get_execution_statistics()

print(f"Executed {stats['total_fragments']} fragments")
print(f"Across {stats['levels_executed']} levels")
print(f"Parallel time: {stats['parallel_time']:.3f}s")
print(f"Estimated speedup: {stats['speedup_factor']:.2f}x")
```

## Validation Methods

### `validate_parallel_correctness(parallel_results: Dict, sequential_results: Dict, tolerance: float = 1e-12) -> bool`

Validate parallel results against sequential results.

**Parameters:**
- `parallel_results` (`Dict[str, AtomicResult]`): Results from parallel execution
- `sequential_results` (`Dict[str, AtomicResult]`): Reference results from sequential execution
- `tolerance` (`float`): Numerical tolerance for comparison (default: 1e-12)

**Returns:**
- `bool`: True if parallel results are within tolerance

**Raises:**
- `ValueError`: If results don't match within tolerance or have structural differences

**Example:**
```python
# Run sequential reference
config_seq = ParallelConfig(execution_mode="serial", max_workers=1)
executor_seq = ParallelManyBodyExecutor(core, config_seq)
sequential_results = executor_seq.execute_full_calculation()

# Run parallel calculation
config_par = ParallelConfig(execution_mode="threading", max_workers=4)
executor_par = ParallelManyBodyExecutor(core, config_par)
parallel_results = executor_par.execute_full_calculation()

# Validate correctness
is_correct = executor_par.validate_parallel_correctness(
    parallel_results, sequential_results, tolerance=1e-12
)

print(f"Parallel correctness: {'✅ PASSED' if is_correct else '❌ FAILED'}")
```

## Properties and Attributes

### Instance Attributes

- `core` (`ManyBodyCore`): The associated ManyBodyCore instance
- `config` (`ParallelConfig`): Parallel execution configuration
- `execution_stats` (`Dict`): Real-time execution statistics

### Private Attributes

- `_dependency_graph` (`NBodyDependencyGraph`): P1-002 dependency graph instance

## Error Handling

### Common Exceptions

**RuntimeError**: Most common exception, raised for:
- Missing P1-002 foundation methods
- Fragment calculation failures
- Parallel execution errors
- Invalid results or timeouts

**ValueError**: Raised for:
- Invalid configuration parameters
- Validation failures
- Structural mismatches in results

### Error Recovery

The executor provides robust error handling:

1. **Fragment Failures**: Individual fragment failures stop level execution
2. **Level Failures**: Level failures stop the entire calculation
3. **Timeout Handling**: Configurable timeouts for long calculations
4. **Resource Management**: Automatic cleanup of parallel resources

### Debugging Tips

1. **Use Serial Mode**: Start with `execution_mode="serial"` for debugging
2. **Check Logs**: Enable debug logging for detailed execution traces
3. **Validate Configuration**: Use small test systems first
4. **Monitor Resources**: Check memory and CPU usage during execution

## Performance Considerations

### Optimal Configuration

- **Threading vs Multiprocessing**: Threading often better for I/O-bound QC calculations
- **Worker Count**: Start with CPU core count, adjust based on memory usage
- **Memory Limits**: Set appropriate per-worker memory limits
- **Timeout Values**: Set reasonable timeouts based on expected calculation time

### Scalability Guidelines

- **Small Systems** (≤3 fragments): Limited parallelization benefit
- **Medium Systems** (4-8 fragments): Good speedup potential
- **Large Systems** (>8 fragments): Excellent parallelization efficiency

### Memory Management

The executor automatically manages memory:
- Per-worker memory limits enforced
- Automatic cleanup between calculations
- QCEngine memory configuration handled
- Memory monitoring in execution statistics

---

For complete examples and usage patterns, see the [Usage Guide](../usage/README.md).