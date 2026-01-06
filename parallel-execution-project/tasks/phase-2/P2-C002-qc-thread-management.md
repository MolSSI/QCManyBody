# QC Thread Management System

**Task ID**: P2-C002
**Task Name**: QC Thread Management System
**Phase**: 2
**Category**: Hybrid Parallelism
**Owner**: Lead Developer
**Estimated Effort**: 3 days
**Priority**: P1 (High)
**Status**: NOT_STARTED

## Description
Implement a comprehensive thread management system that coordinates both fragment-level parallelism and QC-calculation-level parallelism, allowing individual quantum chemistry calculations to use multiple threads while maintaining overall system resource control.

## Acceptance Criteria
- [ ] Configure QC programs to use specified number of threads per calculation
- [ ] Coordinate QC threads with fragment worker threads to prevent oversubscription
- [ ] Support dynamic thread allocation based on available resources
- [ ] Integrate with QCEngine thread configuration
- [ ] Provide thread affinity management when beneficial
- [ ] Monitor and report thread utilization statistics
- [ ] Handle QC program-specific threading requirements (Psi4, NWChem, etc.)
- [ ] Graceful handling of QC programs that don't support threading

## Technical Requirements
- Integration with existing `ParallelManyBodyExecutor`
- Support for major QC programs (Psi4, NWChem, CFOUR, etc.)
- Thread safety for concurrent QC calculations
- Resource conflict detection and prevention
- Platform-specific optimization (Linux, macOS, Windows)
- Memory affinity coordination with thread affinity

## Dependencies
### Prerequisite Tasks
- [x] Phase 1: Core parallel execution infrastructure
- [ ] P2-A001: System Resource Detection Framework
- [ ] P2-C001: Hybrid Parallel Architecture Design

### External Dependencies
- [ ] QCEngine thread configuration capabilities
- [ ] QC program threading documentation and testing
- [ ] Thread affinity libraries (optional, for optimization)

## Deliverables
1. **Primary Deliverable**: `QCThreadManager` class:
   ```python
   class QCThreadManager:
       def configure_qc_threads(self,
                               qc_program: str,
                               threads_per_calc: int,
                               total_workers: int) -> QCThreadConfig:
           """Configure QC program threading parameters."""

       def allocate_thread_resources(self,
                                   fragment_workers: int,
                                   qc_threads_per_worker: int) -> ThreadAllocation:
           """Allocate thread resources across workers."""

       def set_thread_affinity(self, worker_id: int, thread_ids: List[int]):
           """Set CPU affinity for worker threads."""

       def monitor_thread_usage(self) -> ThreadUsageStats:
           """Monitor real-time thread utilization."""
   ```

2. **Configuration Classes**:
   ```python
   @dataclass
   class HybridParallelConfig(ParallelConfig):
       """Extended parallel config with hybrid threading support."""
       qc_threads_per_fragment: int = 1
       enable_thread_affinity: bool = False
       thread_allocation_strategy: str = "balanced"
       max_total_threads: Optional[int] = None

   @dataclass
   class QCThreadConfig:
       """QC program-specific thread configuration."""
       program_name: str
       threads_per_calculation: int
       environment_variables: Dict[str, str]
       qcengine_task_config: Dict[str, Any]

   @dataclass
   class ThreadAllocation:
       """Thread allocation plan for hybrid execution."""
       fragment_workers: int
       qc_threads_per_worker: int
       total_threads_used: int
       cpu_affinity_map: Dict[int, List[int]]
   ```

3. **QC Program Integration**:
   ```python
   class QCProgramThreading:
       """QC program-specific threading implementations."""

       def configure_psi4_threading(self, threads: int) -> Dict[str, Any]:
           """Configure Psi4 threading parameters."""

       def configure_nwchem_threading(self, threads: int) -> Dict[str, Any]:
           """Configure NWChem threading parameters."""

       def configure_cfour_threading(self, threads: int) -> Dict[str, Any]:
           """Configure CFOUR threading parameters."""
   ```

4. **Documentation**: Complete threading guide with QC program examples
5. **Tests**: Comprehensive threading validation and performance tests

## Implementation Approach
### High-Level Steps
1. **QC Program Analysis**: Research threading capabilities of major QC programs
2. **Thread Coordination Design**: Design resource allocation algorithms
3. **QCEngine Integration**: Implement QCEngine thread configuration
4. **Thread Affinity**: Implement CPU affinity management (optional optimization)
5. **Resource Monitoring**: Implement thread usage tracking
6. **QC Program Support**: Implement program-specific threading
7. **Testing and Validation**: Comprehensive testing across QC programs

### Technical Considerations
- QC programs use different threading models (OpenMP, MKL, custom)
- Thread oversubscription can severely hurt performance
- Memory bandwidth competition between threads
- NUMA considerations for thread placement
- QCEngine task configuration integration
- Environment variable management for QC programs

## QC Program Threading Support
### Psi4 Threading
```python
# Psi4 threading configuration
psi4_config = {
    "environment": {
        "OMP_NUM_THREADS": str(qc_threads),
        "MKL_NUM_THREADS": str(qc_threads),
        "PSI_NTHREAD": str(qc_threads)
    },
    "keywords": {
        "num_threads": qc_threads
    }
}
```

### NWChem Threading
```python
# NWChem threading configuration
nwchem_config = {
    "environment": {
        "OMP_NUM_THREADS": str(qc_threads),
        "NWCHEM_NPROC": str(qc_threads)
    },
    "task_config": {
        "ncores": qc_threads
    }
}
```

## Threading Strategies
### Strategy 1: Balanced Threading
- Equal threads per fragment worker
- Total threads = fragment_workers × qc_threads_per_worker
- Simple and predictable resource usage

### Strategy 2: Adaptive Threading
- More threads for larger/complex fragments
- Dynamic allocation based on fragment characteristics
- Optimal resource utilization

### Strategy 3: Memory-Aware Threading
- Consider memory bandwidth limitations
- Reduce threads when memory-bound
- Balance compute and memory resources

## Testing Strategy
### Unit Tests
- [ ] Thread allocation algorithm validation
- [ ] QC program configuration generation
- [ ] Resource conflict detection
- [ ] Thread affinity management
- [ ] Environment variable handling

### Integration Tests
- [ ] QCEngine integration with threading
- [ ] Multi-QC program threading support
- [ ] Fragment worker coordination
- [ ] Resource monitoring accuracy
- [ ] Performance regression testing

### Performance Tests
- [ ] Threading performance on various hardware
- [ ] Comparison: fragment-only vs. hybrid parallelism
- [ ] Thread scaling efficiency measurement
- [ ] Memory bandwidth utilization
- [ ] QC program-specific optimization validation

## Performance Requirements
- **Thread allocation overhead**: <10ms per configuration
- **Resource monitoring overhead**: <1% of total execution time
- **Thread efficiency**: >90% CPU utilization with optimal configuration
- **Memory efficiency**: No significant memory overhead per thread

## Error Handling Strategy
- **QC program detection**: Graceful handling of unsupported programs
- **Thread limit validation**: Prevent oversubscription
- **Affinity failures**: Fallback when affinity setting fails
- **Resource conflicts**: Detection and resolution of resource conflicts

## Integration Points
- **ParallelManyBodyExecutor**: Core execution engine enhancement
- **ParallelConfig**: Extended configuration with hybrid options
- **QCEngine**: Thread configuration for quantum chemistry calculations
- **System Resource Detection**: Optimal thread count determination

## Performance Optimization Opportunities
### CPU Affinity
- Pin QC calculation threads to specific CPU cores
- Reduce cache misses and memory access conflicts
- Optimize for NUMA topology

### Memory Affinity
- Coordinate memory allocation with CPU placement
- Minimize cross-socket memory access
- Optimize for memory bandwidth

### Load Balancing
- Dynamic thread reallocation based on workload
- Consider fragment complexity for thread assignment
- Adaptive scheduling based on real-time performance

## Definition of Done
- [ ] All acceptance criteria met and verified
- [ ] `QCThreadManager` and supporting classes implemented
- [ ] QC program-specific threading support (Psi4, NWChem minimum)
- [ ] Integration with existing parallel execution system
- [ ] Hybrid parallelism showing >20% performance improvement
- [ ] Thread management prevents resource conflicts
- [ ] Unit and integration tests passing (>95% coverage)
- [ ] Performance requirements met
- [ ] Documentation complete with QC program examples
- [ ] Code reviewed and approved
- [ ] No regressions in existing functionality

## Notes & Comments
This task is critical for achieving the performance potential of hybrid parallelism. Key considerations:

**Thread Coordination Complexity**:
- Must carefully balance fragment workers and QC threads
- Avoid oversubscription which can cause severe performance degradation
- Consider memory bandwidth limitations

**QC Program Diversity**:
- Different QC programs use different threading models
- Some programs have optimal thread counts for given problem sizes
- Thread configuration methods vary significantly

**Performance Trade-offs**:
- More QC threads may reduce fragment-level parallelism
- Optimal configuration depends on hardware and calculation characteristics
- Need empirical testing to validate theoretical optimizations

## Timeline
- **Start Date**: Day 5 of Phase 2 (after P2-A001, P2-C001)
- **Target Completion**: Day 8 of Phase 2
- **Critical Path**: Yes (foundational for hybrid parallelism)

## Change Log
| Date | Change | Reason |
|------|--------|--------|
| 2024-09-26 | Task created | Phase 2 planning |