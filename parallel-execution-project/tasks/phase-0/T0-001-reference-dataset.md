# Generate Golden Reference Dataset

**Task ID**: T0-001
**Task Name**: Generate Golden Reference Dataset
**Phase**: 0 (Testing Foundation)
**Owner**: QA Engineer
**Estimated Effort**: 3 days
**Priority**: CRITICAL
**Status**: NOT_STARTED

## Description
Generate a comprehensive golden reference dataset using the current sequential QCManyBody code before any parallel modifications. This dataset will serve as the absolute truth for all regression testing throughout parallel development.

## Acceptance Criteria
- [ ] Reference dataset covers all existing QCManyBody test scenarios
- [ ] Additional test cases added for parallel-specific validation
- [ ] All BSSE types (cp, nocp, vmfc) included
- [ ] All N-body levels (1-4) covered
- [ ] All property types (energy, gradient, hessian) captured
- [ ] Multi-level calculations included
- [ ] Edge cases (supersystem, embedding charges) covered
- [ ] Reference data compressed and stored efficiently
- [ ] Validation script confirms dataset completeness

## Technical Requirements
- Must use UNMODIFIED sequential QCManyBody code (current main branch)
- Generate references for existing test systems: he4, h2o3, ar5, etc.
- Include additional systematic test cases for comprehensive coverage
- Store references in compressed JSON format (zstandard)
- Numerical precision must be at least 1e-14 (working precision)
- Include metadata: QCManyBody version, QC program versions, system specs

## Dependencies
### Prerequisite Tasks
- None (foundational task)

### External Dependencies
- [ ] Access to QC programs (Psi4, NWChem) for reference generation
- [ ] Sufficient computational resources for reference calculations
- [ ] Clean QCManyBody installation without parallel modifications

## Deliverables
1. **Primary Deliverable**: Complete reference dataset in `qcmanybody/tests/reference_data_parallel/`
   - `reference_data_v1.0.json.zst` - Compressed reference dataset
   - `reference_metadata.json` - Dataset metadata and version info
   - `reference_index.json` - Index of all test cases included

2. **Validation Tools**:
   - `scripts/generate_reference_data.py` - Reference generation script
   - `scripts/validate_reference_data.py` - Dataset validation utility
   - `qcmanybody/testing/reference_manager.py` - Reference data access API

3. **Documentation**:
   - Reference dataset documentation with test case descriptions
   - Usage guide for accessing and validating references
   - Versioning strategy for reference data updates

## Implementation Approach
### High-Level Steps
1. **Test Case Inventory**: Catalog all existing QCManyBody tests and identify gaps
2. **Additional Test Design**: Create additional test cases for comprehensive coverage
3. **Reference Generation**: Execute all test cases with sequential code
4. **Data Organization**: Structure and compress reference data efficiently
5. **Validation**: Verify completeness and accuracy of generated references
6. **Documentation**: Document dataset structure and usage

### Test Case Matrix Design
```python
# Systematic test case generation
test_matrix = {
    "systems": ["he2", "he4", "h2o3", "ar5", "mixed_6"],
    "bsse_types": ["cp", "nocp", "vmfc"],
    "max_nbody": [2, 3, 4],
    "drivers": ["energy", "gradient", "hessian"],
    "levels": ["single", "multi"],
    "special_cases": ["supersystem", "embedding_charges"]
}
# Total combinations: 5×3×3×3×2 + special = ~270 test cases
```

### Data Storage Format
```python
# Reference data structure
reference_entry = {
    "test_id": "he4_cp_energy_maxnb3",
    "system": "he4_tetramer",
    "config": {
        "bsse_type": "cp",
        "max_nbody": 3,
        "driver": "energy",
        "levels": {1: "hf", 2: "mp2", 3: "ccsd"}
    },
    "qcmanybody_version": "0.3.1",
    "qc_program": "psi4",
    "qc_version": "1.9.1",
    "timestamp": "2024-09-25T00:00:00Z",
    "results": {
        "return_result": -11.526123456789,
        "properties": {...},
        "component_properties": {...}
    },
    "numerical_precision": 1e-14
}
```

## Testing Strategy
### Reference Generation Validation
- [ ] Compare generated references against existing static test data
- [ ] Verify numerical precision meets requirements
- [ ] Check completeness against test case matrix
- [ ] Validate data integrity after compression/decompression

### Quality Assurance
- [ ] Independent verification of critical test cases
- [ ] Cross-validation with multiple QC programs where possible
- [ ] Statistical analysis of numerical precision
- [ ] Documentation review for completeness

## Definition of Done
- [ ] All acceptance criteria met
- [ ] Reference dataset generated and validated
- [ ] Compression and storage format finalized
- [ ] Access API implemented and tested
- [ ] Documentation complete and reviewed
- [ ] Dataset integrity verified
- [ ] Backup and version control established

## Notes & Comments
This is the **most critical task** for the entire parallel execution project. The quality and completeness of this reference dataset will determine the success of all subsequent regression testing.

**Critical Success Factors:**
- Must use absolutely unmodified sequential code
- Numerical precision must be at working precision (1e-14+)
- Coverage must be comprehensive enough to catch any parallel implementation issues
- Data format must be future-proof and extensible

**Risk Mitigation:**
- Generate references on multiple systems to verify reproducibility
- Create both compressed and uncompressed versions for debugging
- Include checksums and integrity verification
- Document exact software versions and system configurations used

## Timeline
- **Start Date**: TBD (immediately upon project kickoff)
- **Target Completion**: 3 days after start
- **Actual Completion**: TBD

## Change Log
| Date | Change | Reason |
|------|--------|--------|
| 2024-09-25 | Task created | Testing foundation requirement |