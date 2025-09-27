# Fix Threading-based QCEngine Integration

**Task ID**: P1A-001
**Task Name**: Fix Threading-based QCEngine Integration
**Phase**: 1a (Parallelism Fixes)
**Owner**: Lead Developer
**Estimated Effort**: 2 days
**Priority**: HIGH
**Status**: NOT_STARTED

## Description
Resolve the QCEngine global configuration issues that prevent proper execution in multithreaded environments. The current implementation fails with `KeyError: 'ncores'` and `KeyError: 'memory'` when QCEngine attempts to access global configuration values in worker threads.

## Acceptance Criteria
- [ ] Threading-based parallel execution works without QCEngine global config errors
- [ ] All QCEngine task_config parameters are properly initialized in worker threads
- [ ] No dependency on QCEngine global state during parallel fragment execution
- [ ] Backward compatibility maintained with existing QCEngine workflows
- [ ] Performance equivalent to sequential execution per fragment
- [ ] Error handling preserves original QCEngine error messages when appropriate

## Technical Requirements
- Must work with QCEngine v0.33.0+ threading execution model
- Compatible with all QC programs supported by QCEngine (Psi4, NWChem, CFOUR)
- Thread-safe QCEngine configuration management
- No modification of QCEngine library itself (external dependency)
- Preserve existing task_config parameter precedence and overrides

## Dependencies
### Prerequisite Tasks
- [x] P1-001: Dependency graph foundation (completed)
- [x] P1-002: Level-by-level iteration (completed)

### External Dependencies
- [ ] QCEngine v0.33.0+ installed and working
- [ ] At least one QC program (Psi4) available for testing
- [ ] Understanding of QCEngine's internal configuration system

## Deliverables
1. **Primary Deliverable**: Updated `execute_fragment()` method in `qcmanybody/parallel.py`
2. **Configuration Handler**: QCEngine thread-safe initialization logic
3. **Tests**: Unit tests demonstrating threading works with QCEngine
4. **Documentation**: Updated threading execution documentation

## Implementation Approach
### High-Level Steps
1. **Analysis**: Investigate QCEngine's global configuration system and thread safety
2. **Initialization Strategy**: Implement per-thread QCEngine configuration setup
3. **Task Config Enhancement**: Ensure complete task_config parameters for isolation
4. **Error Handling**: Robust fallback when global config unavailable
5. **Testing**: Validate threading with real QC calculations
6. **Documentation**: Update API docs and troubleshooting guides

### Technical Considerations
- QCEngine uses global `_global_values` dict that's not thread-safe
- `get_node_descriptor()` always tries to access globals even with task_config
- Worker threads need complete configuration or proper initialization
- Alternative: Force QCEngine initialization in each worker thread
- Fallback strategy when global config cannot be established

## Testing Strategy
### Unit Tests
- [ ] Threading execution without QCEngine (placeholder mode)
- [ ] Threading with QCEngine using simple HF/STO-3G calculations
- [ ] Multiple concurrent QCEngine calls with different task_config values
- [ ] Error recovery when QCEngine global config fails
- [ ] Configuration parameter precedence (task_config vs globals)

### Integration Tests
- [ ] Full water dimer calculation using threading mode
- [ ] Different QC programs (Psi4 primary, others if available)
- [ ] Various memory and ncores configurations
- [ ] Mixed BSSE types with threading execution
- [ ] Large fragment count systems (stress test)

### Performance Tests
- [ ] Threading overhead compared to sequential execution
- [ ] Memory usage with multiple QCEngine threads
- [ ] Execution time scaling with worker count

## Definition of Done
- [ ] All acceptance criteria met
- [ ] Code reviewed and approved
- [ ] Unit tests written and passing (>95% coverage for new code)
- [ ] Integration tests passing with real QC calculations
- [ ] No threading-related QCEngine errors in test suite
- [ ] Documentation updated in docstrings and user guides
- [ ] Performance regression tests passing

## Notes & Comments
This task addresses the immediate threading compatibility issue discovered during parallel execution testing. The root cause is QCEngine's reliance on global configuration state that isn't properly initialized in worker threads.

Key design decisions:
- Whether to modify QCEngine initialization in workers vs. providing complete task_config
- Error handling strategy when QCEngine global config unavailable
- Thread safety considerations for QCEngine resource management

The solution should be minimal and focused, avoiding complex QCEngine internals modifications.

## Timeline
- **Start Date**: TBD
- **Target Completion**: TBD
- **Actual Completion**: TBD

## Change Log
| Date | Change | Reason |
|------|--------|--------|
| 2024-09-27 | Task created | Threading QCEngine integration issues discovered |