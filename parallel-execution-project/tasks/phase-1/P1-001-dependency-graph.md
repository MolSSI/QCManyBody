# Create N-Body Dependency Graph

**Task ID**: P1-001
**Task Name**: Create N-Body Dependency Graph
**Phase**: 1
**Owner**: Lead Developer
**Estimated Effort**: 5 days
**Priority**: HIGH
**Status**: NOT_STARTED

## Description
Implement the core dependency management system that can extract N-body levels from fragment labels and create a dependency graph ensuring proper execution order (monomers → dimers → trimers → N-mers).

## Acceptance Criteria
- [ ] Function `extract_nbody_level(label: str) -> int` correctly identifies N-body level from any fragment label
- [ ] Function `group_by_level(fragments) -> Dict[int, List[Fragment]]` groups fragments by dependency level
- [ ] Dependency validation ensures all required lower-level fragments are present
- [ ] Works with all BSSE treatment types (cp, nocp, vmfc)
- [ ] Handles multi-level calculations correctly
- [ ] Full unit test coverage (>95%)

## Technical Requirements
- Must parse existing QCManyBody fragment labels (format: `["(method)", [real_atoms], [basis_atoms]]`)
- Compatible with existing `labeler()` and `delabeler()` utilities in `qcmanybody/utils.py`
- Thread-safe for parallel execution
- Memory efficient for large fragment sets (100+ fragments)
- Integration with existing `ManyBodyCore.compute_map` structure

## Dependencies
### Prerequisite Tasks
- None (foundational task)

### External Dependencies
- [ ] Understanding of existing fragment labeling system
- [ ] Access to `qcmanybody/utils.py` labeler functions
- [ ] Example fragment labels from existing test cases

## Deliverables
1. **Primary Deliverable**: `qcmanybody/dependency.py` module with core classes:
   - `NBodyDependencyGraph`
   - `FragmentDependency`
   - Utility functions for level extraction and grouping
2. **Documentation**: API documentation with examples
3. **Tests**: Comprehensive unit tests in `qcmanybody/tests/test_dependency.py`
4. **Examples**: Usage examples for different fragment types

## Implementation Approach
### High-Level Steps
1. **Analysis**: Study existing fragment labeling in `core.py` and `utils.py`
2. **Design**: Create dependency graph data structure
3. **Implementation**: Core dependency management functions
4. **Integration**: Integrate with existing `iterate_molecules()` method
5. **Testing**: Comprehensive test suite
6. **Documentation**: API docs and usage examples

### Technical Considerations
- Fragment labels use format `["method", [real_atoms], [basis_atoms]]`
- N-body level = length of `real_atoms` list
- Must handle ghost atoms in `basis_atoms` correctly
- Consider edge cases: supersystem calculations, embedding charges
- Performance optimization for large dependency graphs

## Testing Strategy
### Unit Tests
- [ ] Extract N-body level from various fragment label formats
- [ ] Group fragments correctly by dependency level
- [ ] Validate dependency chains for different BSSE types
- [ ] Handle edge cases (empty fragments, supersystem)
- [ ] Performance test with large fragment sets

### Integration Tests
- [ ] Integration with existing `ManyBodyCore.iterate_molecules()`
- [ ] Compatibility with all BSSE treatment types
- [ ] Multi-level calculation scenarios
- [ ] Real test case data from existing test suite

### Performance Tests
- [ ] Memory usage with 100+ fragments
- [ ] Execution time for dependency graph construction
- [ ] Scalability analysis

## Definition of Done
- [ ] All acceptance criteria met
- [ ] Code reviewed and approved
- [ ] Unit tests written and passing (>95% coverage)
- [ ] Integration tests passing
- [ ] Documentation updated in docstrings
- [ ] Performance benchmarks completed
- [ ] Compatible with existing QCManyBody test suite

## Notes & Comments
This is the foundational task for parallel execution. The dependency graph system must be robust and correct as all subsequent parallel execution depends on it.

Key design decisions:
- Whether to modify existing `iterate_molecules()` or create new interface
- How to handle supersystem calculations in dependency graph
- Performance vs. memory trade-offs for large systems

## Timeline
- **Start Date**: TBD
- **Target Completion**: TBD
- **Actual Completion**: TBD

## Change Log
| Date | Change | Reason |
|------|--------|--------|
| 2024-09-24 | Task created | Initial planning |