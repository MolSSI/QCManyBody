# Known Issues - QCManyBody Parallel Execution

## Issue #001: Analysis Phase Data Structure Mismatch

**Status**: IDENTIFIED - Needs Resolution
**Priority**: MEDIUM
**Discovered**: 2024-09-27 during P1A-001 testing

### Description
After successful parallel execution of fragment calculations, the analysis phase fails in `core.analyze()` with:
```
AttributeError: 'AtomicResult' object has no attribute 'keys'
File "qcmanybody/core.py", line 643, in analyze
available_properties.update(property_data.keys())
```

### Root Cause Analysis
- **Location**: `qcmanybody/core.py:643` in `analyze()` method
- **Issue**: Code expects `property_data` to be a dictionary with `.keys()` method
- **Reality**: `property_data` is an `AtomicResult` object, not a dictionary
- **Context**: This occurs after successful parallel fragment execution

### Impact
- **Severity**: Does not affect parallel execution correctness
- **Scope**: Prevents final analysis and result aggregation
- **Workaround**: Parallel execution works; analysis structure needs updating

### Technical Details
```python
# Failing code in core.py:643
available_properties.update(property_data.keys())  # AtomicResult has no .keys()

# Expected: property_data should be Dict[str, Any]
# Actual: property_data is AtomicResult object
```

### Resolution Plan
- **Phase**: Post P1A-002 (not blocking multiprocessing development)
- **Approach**: Update analysis phase to handle AtomicResult objects correctly
- **Files**: `qcmanybody/core.py` analysis methods
- **Testing**: Validate with both sequential and parallel execution results

### Reproduction
1. Run water4 cluster test with threading mode
2. Parallel execution completes successfully
3. Analysis phase fails with AttributeError

---

**Next Issue ID**: #002