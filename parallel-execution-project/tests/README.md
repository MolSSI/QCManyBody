# Testing Strategy for Parallel Execution Development

## Overview

This directory contains the comprehensive testing strategy for implementing parallel execution in QCManyBody while ensuring **exact numerical reproduction** of the original sequential code.

## 🎯 Core Testing Principles

1. **Zero Tolerance for Numerical Differences**: Parallel results must be identical to sequential results within 1e-12 tolerance
2. **Continuous Validation**: Every code change must pass regression tests
3. **Multi-Dimensional Testing**: Test across system sizes, BSSE types, worker counts, and QC programs
4. **Fail-Fast Strategy**: Catch regressions immediately during development

## 📁 Testing Documentation

### Core Documents
- **[regression-testing-strategy.md](regression-testing-strategy.md)**: Comprehensive strategy for validating parallel vs sequential results
- **[ci-testing-approach.md](ci-testing-approach.md)**: Continuous integration and automated testing approach

## 🧪 Testing Categories

### 1. Golden Reference Tests (Foundational)
**Purpose**: Establish numerical baselines before any parallel modifications

```python
# Generate comprehensive reference dataset
python scripts/generate_reference_data.py

# Includes:
- All existing QCManyBody test cases
- Additional parallel-specific test scenarios
- Multi-dimensional test matrix (BSSE × N-body × system size)
```

**Critical**: Must be generated with unmodified sequential code

### 2. Regression Tests (Development)
**Purpose**: Validate every parallel modification reproduces original results

```python
# Run after every code change
pytest qcmanybody/tests/test_parallel_regression.py -v

# Tests:
- Dependency graph construction preserves iteration order
- Parallel execution matches sequential (all worker counts)
- Load balancing strategies produce identical results
```

**Tolerance**: 1e-12 absolute and relative tolerance

### 3. Integration Tests (QC Programs)
**Purpose**: Ensure compatibility with all supported quantum chemistry programs

```python
# Test with real QC calculations
pytest qcmanybody/tests/test_computer_*.py --parallel-mode=multiprocessing
pytest qcmanybody/tests/test_examples.py --parallel-mode=multiprocessing

# Validates:
- Psi4, NWChem, CFOUR compatibility
- Process isolation and resource management
- Error handling consistency
```

### 4. Performance Tests (Validation)
**Purpose**: Verify parallel execution improves performance while maintaining correctness

```python
# Performance benchmarks with correctness validation
python scripts/benchmark_parallel_execution.py --validate-accuracy

# Measures:
- Speedup factors (2-6× target)
- Memory usage efficiency
- Load balancing effectiveness
```

## 🔧 Testing Tools and Utilities

### Validation Framework
```python
from qcmanybody.testing import ParallelRegressionTester

tester = ParallelRegressionTester(tolerance=1e-12)
assert tester.validate_parallel_result(
    parallel_result,
    reference_key="he4_cp_energy"
)
```

### Reference Data Management
```python
from qcmanybody.testing import ReferenceDataManager

manager = ReferenceDataManager()
reference = manager.load_reference("h2o3_vmfc_gradient")
manager.validate_against_reference(test_result, reference)
```

### CI Integration
```python
# Automated CI testing with GitHub Actions
# Multi-stage pipeline:
# 1. Quick validation (< 5 min) - unit tests, dependency tests
# 2. Parallel validation (< 15 min) - regression tests
# 3. QC integration (< 30 min) - real calculations
# 4. Performance regression (main branch only)
```

## 📊 Testing Matrix

### Test Dimensions
| Dimension | Values | Purpose |
|-----------|--------|---------|
| **System Size** | 2, 3, 4, 6, 8 fragments | Scalability |
| **BSSE Type** | cp, nocp, vmfc | All correction methods |
| **N-body Level** | 2, 3, 4 | Dependency chains |
| **Workers** | 1, 2, 4, 8 | Parallel scaling |
| **QC Program** | Psi4, NWChem, CFOUR | Program compatibility |
| **Driver** | energy, gradient, hessian | Property types |

### Example Test Case Matrix
```python
@pytest.mark.parametrize("fragments", [2, 3, 4, 6])
@pytest.mark.parametrize("bsse_type", ["cp", "nocp", "vmfc"])
@pytest.mark.parametrize("workers", [1, 2, 4])
def test_comprehensive_regression_matrix(fragments, bsse_type, workers):
    """Test all combinations for regression."""
    # 2×3×3 = 18 test combinations per system
```

## 🚀 Development Workflow

### Before Implementation
1. **Generate References**: Create golden dataset with current sequential code
2. **Baseline Tests**: Ensure all existing tests pass
3. **Environment Setup**: Configure development and CI environments

### During Development
1. **Unit Tests**: Test individual components (dependency graph, load balancing)
2. **Integration Tests**: Test parallel execution with simple cases
3. **Regression Tests**: Validate against reference data continuously

### Before Merge
1. **Full Regression Suite**: All reference tests must pass
2. **Multi-Worker Validation**: Test with 1, 2, 4 workers
3. **QC Program Testing**: At least one QC program must work
4. **Performance Check**: No significant performance regression

### Production Ready
1. **Comprehensive Testing**: All QC programs, all test scenarios
2. **Stress Testing**: Large systems, edge cases
3. **Long-term Stability**: Extended test runs
4. **Documentation**: User guides and examples

## ⚠️ Critical Testing Requirements

### Mathematical Correctness
- **Numerical Identity**: Results must be identical within working precision
- **BSSE Treatment**: All basis set superposition error corrections must work
- **Multi-level**: Complex level calculations must be exact
- **Edge Cases**: Supersystem, embedding charges, etc.

### Error Handling
- **Consistent Failures**: Parallel and sequential must fail the same way
- **Resource Management**: Proper cleanup on success and failure
- **Timeout Handling**: Graceful handling of stuck calculations

### Performance Requirements
- **Speedup Validation**: Must achieve target speedups (2-6×)
- **Memory Efficiency**: Reasonable memory usage per worker
- **Scalability**: Performance should scale with worker count

## 🔍 Debugging and Troubleshooting

### Common Issues and Solutions

#### Numerical Differences
```python
# Debugging numerical mismatches
def debug_numerical_difference(par_result, seq_result):
    diff = abs(par_result - seq_result)
    if diff > 1e-12:
        print(f"NUMERICAL MISMATCH: {par_result} vs {seq_result}")
        print(f"Absolute difference: {diff:.2e}")
        # Check:
        # - Fragment ordering consistency
        # - Floating point precision issues
        # - QC program determinism
```

#### Worker Process Issues
```python
# Debugging worker failures
def debug_worker_failure(worker_logs):
    # Check:
    # - Process initialization
    # - Resource conflicts
    # - QC program compatibility
    # - Memory limitations
```

#### Load Balancing Problems
```python
# Debugging load imbalance
def debug_load_balancing(worker_times):
    utilization = np.std(worker_times) / np.mean(worker_times)
    if utilization > 0.2:  # >20% variance
        print("LOAD IMBALANCE DETECTED")
        # Check:
        # - Cost estimation accuracy
        # - Work distribution algorithm
        # - System heterogeneity
```

## 📈 Success Metrics

### Correctness Metrics
- [ ] **100%** of reference tests pass
- [ ] **Zero** numerical differences > 1e-12
- [ ] **100%** of existing tests continue to pass
- [ ] **Identical** error handling between parallel/sequential

### Performance Metrics
- [ ] **>2×** speedup on 4-core systems
- [ ] **>4×** speedup on 8-core systems
- [ ] **<20%** memory overhead vs sequential
- [ ] **<15%** worker utilization variance

### Quality Metrics
- [ ] **>95%** test coverage for parallel code
- [ ] **<1%** test failure rate in CI
- [ ] **100%** QC program compatibility
- [ ] **Zero** production regression bugs

## 📞 Support and Resources

### Documentation
- QCManyBody testing documentation: `qcmanybody/tests/README.md`
- Parallel execution API: `parallel-execution-project/docs/api-design.md`
- Architecture overview: `parallel-execution-project/docs/architecture.md`

### Development Tools
- Reference data generator: `scripts/generate_reference_data.py`
- Parallel validation utility: `scripts/validate_parallel_accuracy.py`
- Performance benchmarks: `scripts/benchmark_parallel_execution.py`
- CI failure analyzer: `scripts/analyze_ci_failures.py`

### Contact
- Technical issues: GitHub Issues with `parallel-execution` label
- Test failures: Include full test output and system information
- Performance questions: Include benchmark results and system specs

---

This comprehensive testing strategy ensures that parallel execution development maintains the mathematical rigor essential for quantum chemistry while providing the performance improvements needed for large many-body calculations.