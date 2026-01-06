# Parallel Execution Architecture

## Overview

This document describes the technical architecture for implementing parallel N-body calculations in QCManyBody while maintaining mathematical correctness and backward compatibility.

## Current Architecture Analysis

### Sequential Execution Flow
```
ManyBodyComputer.compute()
└── iterate_molecules() [SEQUENTIAL]
    ├── Fragment 1 (monomer)
    ├── Fragment 2 (monomer)
    ├── Fragment 3 (dimer)
    ├── Fragment 4 (dimer)
    ├── Fragment 5 (trimer)
    └── ...
└── ManyBodyCore.analyze(all_results)
```

### Identified Issues
1. **No dependency ordering**: `iterate_molecules()` provides fragments in arbitrary order
2. **Blocking execution**: Each `qcng.compute()` call blocks until completion
3. **Resource underutilization**: Only one CPU core active during QC calculations
4. **No parallelization interface**: No API for controlling parallel execution

## Proposed Parallel Architecture

### Level-by-Level Execution Flow
```
ParallelManyBodyComputer.compute()
└── iterate_molecules_by_level() [NEW]
    ├── Level 1: [monomer1, monomer2, monomer3] ── PARALLEL ──→ Results1
    ├── Level 2: [dimer1-2, dimer1-3, dimer2-3] ── PARALLEL ──→ Results2
    ├── Level 3: [trimer1-2-3] ─────────────────── PARALLEL ──→ Results3
    └── ...
└── ManyBodyCore.analyze(all_results)
```

### Key Architectural Components

#### 1. Dependency Management
```python
class NBodyDependencyGraph:
    def extract_nbody_level(self, label: str) -> int:
        """Extract N-body level from fragment label"""

    def group_by_level(self, fragments: List[Fragment]) -> Dict[int, List[Fragment]]:
        """Group fragments by N-body dependency level"""

    def validate_dependencies(self, results: Dict[str, Result]) -> bool:
        """Ensure all required lower-level results are available"""
```

#### 2. Parallel Execution Engine
```python
class ParallelExecutionEngine:
    def __init__(self, mode: ParallelMode, max_workers: int):
        self.mode = mode  # MULTIPROCESSING, MPI, THREAD
        self.max_workers = max_workers

    def execute_level(self, fragments: List[Fragment]) -> Dict[str, Result]:
        """Execute all fragments at a given level in parallel"""

    def aggregate_results(self, results: List[Future]) -> Dict[str, Result]:
        """Collect and validate parallel execution results"""
```

#### 3. Load Balancing Strategy
```python
class LoadBalancer:
    def estimate_cost(self, fragment: Fragment) -> float:
        """Estimate computational cost based on fragment size and method"""

    def distribute_work(self, fragments: List[Fragment], workers: int) -> List[List[Fragment]]:
        """Optimally distribute fragments across workers"""
```

## Implementation Strategy

### Phase 1: Foundation
- Implement `NBodyDependencyGraph` for fragment ordering
- Create `iterate_molecules_by_level()` method
- Add basic parallel execution framework

### Phase 2: Core Parallel Engine
- Multiprocessing implementation using `ProcessPoolExecutor`
- Result aggregation and error handling
- Load balancing algorithms

### Phase 3: Advanced Features
- MPI support for HPC clusters
- Adaptive parallelization based on system resources
- Memory optimization for large calculations

## API Design

### New Parallel Interface
```python
from qcmanybody import ParallelManyBodyComputer, ParallelMode

# Basic parallel execution
computer = ParallelManyBodyComputer(
    input_data=mb_input,
    parallel_mode=ParallelMode.MULTIPROCESSING,
    max_workers=4
)
result = computer.compute()

# Advanced configuration
computer = ParallelManyBodyComputer(
    input_data=mb_input,
    parallel_config={
        'mode': ParallelMode.MPI,
        'max_workers': 16,
        'load_balancing': 'dynamic',
        'memory_limit': '8GB'
    }
)
```

### Backward Compatibility
```python
# Existing interface continues to work unchanged
computer = ManyBodyComputer(input_data=mb_input)
result = computer.compute()  # Sequential execution as before
```

## Performance Considerations

### Expected Speedup
- **4-fragment system**: 2-3× speedup on 4-core systems
- **6-fragment system**: 3-5× speedup on 8-core systems
- **8+ fragment systems**: 4-6× speedup on 16+ core systems

### Memory Usage
- Each worker process: ~200-500MB base + QC program memory
- Shared fragment data: Minimal overhead with proper serialization
- Peak memory: Max workers × QC program memory requirements

### Bottlenecks & Mitigation
1. **QC program initialization**: Pre-warm worker processes
2. **Result serialization**: Optimize data transfer protocols
3. **Load imbalance**: Dynamic work stealing algorithms
4. **Memory pressure**: Adaptive worker count based on available RAM

## Error Handling & Reliability

### Failure Recovery
```python
class ParallelExecutionResult:
    successful_fragments: Dict[str, Result]
    failed_fragments: Dict[str, Exception]
    retry_count: int

    def retry_failures(self) -> 'ParallelExecutionResult':
        """Retry failed fragments with different configuration"""
```

### Validation Strategy
- **Result verification**: Compare parallel vs sequential for test cases
- **Checkpointing**: Save intermediate results for long calculations
- **Graceful degradation**: Fall back to sequential on persistent failures

## Integration Points

### QCEngine Compatibility
- Leverage existing `qcng.compute()` interface
- Respect QC program threading limitations
- Handle program-specific parallel configurations

### ManyBodyCore Integration
- Minimal changes to existing analysis code
- Preserve all BSSE treatment options (cp, nocp, vmfc)
- Maintain multi-level calculation support

## Testing Strategy

### Unit Tests
- Dependency graph construction and validation
- Fragment grouping and ordering algorithms
- Load balancing correctness

### Integration Tests
- End-to-end parallel execution with all supported QC programs
- Memory usage and performance benchmarks
- Error handling and recovery scenarios

### Performance Tests
- Scalability analysis across different system configurations
- Comparison with sequential execution for correctness
- Real-world many-body system benchmarks