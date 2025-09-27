# Create N-Body Dependency Graph (Enhanced with Testing)

**Task ID**: P1-001
**Task Name**: Create N-Body Dependency Graph
**Phase**: 1 (Dependency Analysis & Architecture)
**Owner**: Lead Developer
**Estimated Effort**: 5 days
**Priority**: HIGH
**Status**: NOT_STARTED

## Description
Implement the core dependency management system that can extract N-body levels from fragment labels and create a dependency graph ensuring proper execution order (monomers → dimers → trimers → N-mers). This enhanced version includes comprehensive testing requirements to ensure no disruption to existing fragment iteration.

## Acceptance Criteria
- [ ] Function `extract_nbody_level(label: str) -> int` correctly identifies N-body level from any fragment label
- [ ] Function `group_by_level(fragments) -> Dict[int, List[Fragment]]` groups fragments by dependency level
- [ ] Dependency validation ensures all required lower-level fragments are present
- [ ] Works with all BSSE treatment types (cp, nocp, vmfc)
- [ ] Handles multi-level calculations correctly
- [ ] **TESTING**: Fragment ordering preservation validated (identical to original `iterate_molecules()`)
- [ ] **TESTING**: No duplicate molecules across levels
- [ ] **TESTING**: All existing test cases continue to pass
- [ ] Full unit test coverage (>95%)

## Technical Requirements
- Must parse existing QCManyBody fragment labels (format: `["(method)", [real_atoms], [basis_atoms]]`)
- Compatible with existing `labeler()` and `delabeler()` utilities in `qcmanybody/utils.py`
- Thread-safe for parallel execution
- Memory efficient for large fragment sets (100+ fragments)
- Integration with existing `ManyBodyCore.compute_map` structure
- **CRITICAL**: Must not change the set of molecules generated, only their ordering

## Dependencies
### Prerequisite Tasks
- [x] T0-001: Generate Golden Reference Dataset
- [x] T0-002: Implement Reference Validation Framework

### External Dependencies
- [ ] Understanding of existing fragment labeling system
- [ ] Access to `qcmanybody/utils.py` labeler functions
- [ ] Example fragment labels from existing test cases
- [ ] Reference validation framework from T0-002

## Deliverables
1. **Primary Deliverable**: `qcmanybody/dependency.py` module with core classes:
   - `NBodyDependencyGraph` - Main dependency management class
   - `FragmentDependency` - Individual fragment dependency representation
   - Utility functions for level extraction and grouping

2. **Testing Deliverables**:
   - `qcmanybody/tests/test_dependency_graph.py` - Comprehensive unit tests
   - `qcmanybody/tests/test_fragment_ordering.py` - Fragment iteration preservation tests
   - Integration tests with existing test suite

3. **Validation Tools**:
   - `scripts/validate_dependency_graph.py` - Standalone validation utility
   - `scripts/compare_fragment_iteration.py` - Before/after iteration comparison

4. **Documentation**:
   - API documentation with examples
   - Dependency graph theory and implementation notes

## Implementation Approach
### High-Level Steps
1. **Analysis**: Study existing fragment labeling in `core.py` and `utils.py`
2. **Design**: Create dependency graph data structure
3. **Implementation**: Core dependency management functions
4. **Testing**: Comprehensive test suite with fragment ordering validation
5. **Integration**: Integrate with existing `iterate_molecules()` method
6. **Validation**: Regression testing against reference dataset
7. **Documentation**: API docs and usage examples

### Fragment Ordering Validation Strategy
```python
class FragmentOrderingValidator:
    """Validates that dependency graph preserves existing fragment iteration."""

    def __init__(self, reference_tester):
        self.reference_tester = reference_tester

    def validate_fragment_preservation(self, original_core, new_core):
        """Ensure new iteration produces identical fragment set."""

        # Get original fragment set
        original_fragments = {
            (label, mol.symbols, mol.geometry)
            for chem, label, mol in original_core.iterate_molecules()
        }

        # Get new level-by-level fragment set
        new_fragments = set()
        for level, level_fragments in new_core.iterate_molecules_by_level():
            for chem, label, mol in level_fragments:
                new_fragments.add((label, mol.symbols, mol.geometry))

        # Validate identical sets
        if original_fragments != new_fragments:
            missing_in_new = original_fragments - new_fragments
            extra_in_new = new_fragments - original_fragments

            raise ValidationError(
                f"Fragment set mismatch:\n"
                f"Missing in new: {missing_in_new}\n"
                f"Extra in new: {extra_in_new}"
            )

        return True

    def validate_dependency_ordering(self, new_core):
        """Ensure dependency ordering is mathematically correct."""
        for level, level_fragments in new_core.iterate_molecules_by_level():
            for chem, label, mol in level_fragments:
                _, real_atoms, basis_atoms = delabeler(label)
                actual_nbody = len(real_atoms)

                if actual_nbody != level:
                    raise ValidationError(
                        f"Fragment {label} has {actual_nbody}-body but in level {level}"
                    )

        return True
```

### Testing Requirements Enhancement
```python
# Enhanced testing requirements for P1-001

class TestDependencyGraph:
    """Comprehensive dependency graph testing."""

    def test_fragment_level_extraction(self):
        """Test N-body level extraction from all label types."""
        test_cases = [
            ('["hf", [1], [1]]', 1),           # Monomer
            ('["mp2", [1,2], [1,2]]', 2),      # Dimer
            ('["ccsd", [1,2,3], [1,2,3,4]]', 3), # Trimer with ghost
        ]

        for label, expected_level in test_cases:
            assert extract_nbody_level(label) == expected_level

    def test_fragment_grouping_by_level(self):
        """Test fragment grouping preserves all fragments."""
        # Test with real fragment sets from existing tests
        for test_system in ["he4", "h2o3", "ar5"]:
            original_fragments = generate_test_fragments(test_system)
            grouped = group_by_level(original_fragments)

            # Verify all fragments preserved
            total_grouped = sum(len(level_frags) for level_frags in grouped.values())
            assert total_grouped == len(original_fragments)

    def test_regression_against_reference(self):
        """CRITICAL: Test against golden reference dataset."""
        for reference_key in self.reference_tester.get_all_reference_keys():
            if not reference_key.startswith("dependency_"):
                continue

            # Load reference calculation
            ref_result = self.reference_tester.load_reference(reference_key)

            # Perform same calculation with new dependency graph
            new_result = compute_with_dependency_graph(reference_key)

            # Must be identical
            assert self.reference_tester.validate_result(new_result, reference_key)
```

## Enhanced Testing Strategy
### Unit Tests (Comprehensive)
- [ ] Extract N-body level from various fragment label formats
- [ ] Group fragments correctly by dependency level
- [ ] Validate dependency chains for different BSSE types
- [ ] Handle edge cases (empty fragments, supersystem)
- [ ] Performance test with large fragment sets (100+ fragments)

### Regression Tests (Critical)
- [ ] **Fragment Set Preservation**: New iteration produces identical molecule set
- [ ] **Dependency Order Correctness**: All N-body dependencies respected
- [ ] **Reference Dataset Validation**: All reference calculations reproduce exactly
- [ ] **Existing Test Compatibility**: All existing tests continue to pass
- [ ] **Multi-dimensional Testing**: Test across BSSE types, system sizes, N-body levels

### Integration Tests (Essential)
- [ ] Integration with existing `ManyBodyCore.iterate_molecules()`
- [ ] Compatibility with all BSSE treatment types
- [ ] Multi-level calculation scenarios
- [ ] Real test case data from existing test suite

### Performance Tests (Important)
- [ ] Memory usage with 100+ fragments
- [ ] Execution time for dependency graph construction
- [ ] Scalability analysis for large systems

## Definition of Done (Enhanced)
- [ ] All acceptance criteria met
- [ ] Code reviewed and approved
- [ ] Unit tests written and passing (>95% coverage)
- [ ] **Fragment ordering preservation validated**
- [ ] **All regression tests pass with 1e-12 tolerance**
- [ ] **All existing QCManyBody tests continue to pass**
- [ ] Integration tests passing
- [ ] Documentation updated in docstrings
- [ ] Performance benchmarks completed
- [ ] CI/CD pipeline passing

## Notes & Comments
This is the **foundational task** for parallel execution, enhanced with comprehensive testing requirements. The dependency graph system must be robust and correct as all subsequent parallel execution depends on it.

**CRITICAL SUCCESS FACTORS:**
- **Mathematical Correctness**: Dependency ordering must be perfect
- **Fragment Preservation**: Cannot lose or duplicate any molecules
- **Numerical Identity**: All calculations must produce identical results
- **Backward Compatibility**: Existing code must continue to work

**Enhanced Testing Approach:**
- Every change validated against golden reference dataset
- Fragment ordering preservation explicitly tested
- Multi-dimensional regression testing across all combinations
- Performance impact measured and bounded

**Key Design Decisions:**
- Whether to modify existing `iterate_molecules()` or create new interface
- How to handle supersystem calculations in dependency graph
- Performance vs. memory trade-offs for large systems
- Integration strategy with existing compute_map structure

## Timeline
- **Start Date**: TBD (after Phase 0 completion)
- **Target Completion**: 5 days after start (includes comprehensive testing)
- **Actual Completion**: TBD

## Change Log
| Date | Change | Reason |
|------|--------|--------|
| 2024-09-25 | Task enhanced with testing requirements | Integration with testing foundation |
| 2024-09-24 | Task created | Initial planning |