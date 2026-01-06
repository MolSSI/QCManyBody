# Implement Ordered Fragment Iteration

**Task ID**: P1-002
**Task Name**: Refactor iterate_molecules() for Ordered Execution
**Phase**: 1
**Owner**: Core Developer
**Estimated Effort**: 3 days
**Priority**: HIGH
**Status**: ✅ COMPLETED (2024-09-26)

## Description
Create a new method `iterate_molecules_by_level()` that provides fragments grouped by N-body dependency level, enabling level-by-level parallel execution while maintaining backward compatibility.

## Acceptance Criteria
- [x] New method `iterate_molecules_by_level()` returns `Iterator[Tuple[int, str, str, Molecule]]`
- [x] Groups fragments by N-body level (1-body, 2-body, etc.)
- [x] Maintains exact same fragment generation as existing `iterate_molecules()`
- [x] Preserves backward compatibility - existing code unaffected
- [x] Handles all BSSE types and multi-level calculations
- [x] Eliminates duplicate molecules within and across levels

## Technical Requirements
- Must use dependency graph from P1-001
- Integration with existing `ManyBodyCore` class in `qcmanybody/core.py`
- No changes to existing `iterate_molecules()` method signature or behavior
- Memory efficient - doesn't precompute all molecules at once
- Generator-based implementation for large calculations

## Dependencies
### Prerequisite Tasks
- [x] P1-001: Create N-body dependency graph implementation

### External Dependencies
- [x] Access to existing `ManyBodyCore` class structure
- [x] Understanding of `compute_map` data structure
- [x] Integration with fragment deduplication logic

## Deliverables
1. **Primary Deliverable**: Enhanced `ManyBodyCore` class with new method:
   ```python
   def iterate_molecules_by_level(self) -> Iterator[Tuple[int, List[Tuple[str, str, Molecule]]]]:
       """Iterate molecules grouped by N-body dependency level."""
   ```
2. **Documentation**: Method documentation and usage examples
3. **Tests**: Unit and integration tests
4. **Examples**: Code examples showing level-by-level iteration

## Implementation Approach
### High-Level Steps
1. **Integration**: Import and use dependency graph from P1-001
2. **Grouping Logic**: Group fragments from `compute_map` by N-body level
3. **Deduplication**: Ensure no duplicate molecules within or across levels
4. **Generator**: Implement as memory-efficient generator
5. **Testing**: Validate against existing `iterate_molecules()` output
6. **Documentation**: Add comprehensive docstring and examples

### Technical Considerations
- Use existing `done_molecules` set pattern from current implementation
- Maintain exact same molecule generation logic (fragment extraction, ghost atoms, etc.)
- Consider memory usage for large numbers of fragments
- Handle edge cases: supersystem, embedding charges, multi-level calculations
- Preserve all molecule metadata and properties

## Testing Strategy
### Unit Tests
- [ ] Verify fragments grouped correctly by N-body level
- [ ] Check no duplicate molecules within levels
- [ ] Ensure all fragments from `iterate_molecules()` are included
- [ ] Test with different BSSE types (cp, nocp, vmfc)
- [ ] Edge cases: single fragment, supersystem calculations

### Integration Tests
- [ ] Compare total fragments with existing `iterate_molecules()`
- [ ] Validate molecule properties are identical
- [ ] Test with real QCManyBody calculation inputs
- [ ] Multi-level calculation scenarios
- [ ] Large fragment system performance

### Regression Tests
- [ ] All existing tests continue to pass
- [ ] No behavioral changes to existing functionality
- [ ] Performance impact assessment

## Definition of Done
- [ ] All acceptance criteria met
- [ ] Code reviewed and approved
- [ ] New method implemented and tested
- [ ] Backward compatibility verified
- [ ] Documentation updated
- [ ] Integration tests passing
- [ ] No regression in existing functionality

## Notes & Comments
This task is critical for enabling parallel execution while maintaining compatibility. The implementation must be careful to preserve exact molecule generation behavior.

Key considerations:
- Whether to implement as method of `ManyBodyCore` or separate utility
- How to handle the complex `compute_map` structure efficiently
- Memory vs. performance trade-offs for large calculations
- Integration strategy with parallel execution engine

## Timeline
- **Start Date**: TBD (after P1-001 completion)
- **Target Completion**: TBD
- **Actual Completion**: TBD

## Change Log
| Date | Change | Reason |
|------|--------|--------|
| 2024-09-24 | Task created | Initial planning |