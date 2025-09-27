# Regression Testing Strategy for Parallel Development

## Objective

Ensure that parallel execution modifications reproduce **exactly** the same numerical results as the original sequential code throughout development, preventing regression bugs and maintaining mathematical correctness.

## Current QCManyBody Test Architecture Analysis

### Existing Test Categories
1. **Unit Tests**: Fast, no external deps (`test_utils.py`, `test_schema_keywords.py`)
2. **Static-Data Regression Tests**: Pre-stored results, no QC calculations (`test_core_*.py`)
3. **End-to-End Tests**: Full QC calculations (`test_computer_*.py`)
4. **Examples**: User-friendly demonstrations (`test_examples.py`)

### Key Validation Functions
- `compare_results()` in `qcmanybody/tests/utils.py`: Main result comparison
- `compare()`: Numerical comparison with 1e-7 relative/absolute tolerance
- `qcelemental.testing.compare_values()`: QCElemental's validation utilities

## Comprehensive Regression Testing Strategy

### 1. Golden Reference Generation

Create comprehensive reference datasets before any parallel modifications:

```python
# Generate reference results with current sequential code
class ReferenceGenerator:
    def generate_reference_suite(self):
        """Generate complete reference dataset for regression testing."""

        test_systems = [
            # Small systems (fast validation)
            ("he_dimer", {"max_nbody": 2, "bsse_type": ["cp", "nocp"]}),
            ("h2o_trimer", {"max_nbody": 3, "bsse_type": ["cp", "nocp", "vmfc"]}),

            # Medium systems (thorough testing)
            ("he_tetramer", {"max_nbody": 4, "bsse_type": ["cp", "nocp", "vmfc"]}),
            ("ar_pentamer", {"max_nbody": 3, "bsse_type": ["cp", "vmfc"]}),

            # Multi-level calculations
            ("mixed_levels", {"levels": {1: "hf", 2: "mp2", 3: "ccsd"}}),

            # Edge cases
            ("supersystem_only", {"supersystem_ie_only": True}),
            ("embedding_charges", {"embedding_charges": {...}}),
        ]

        for system_name, config in test_systems:
            for driver in ["energy", "gradient", "hessian"]:
                reference_data = self._compute_reference(system_name, driver, config)
                self._save_reference(f"{system_name}_{driver}", reference_data)
```

### 2. Parallel vs Sequential Validation Framework

```python
class ParallelRegressionTester:
    """Validate parallel execution against sequential reference."""

    def __init__(self, tolerance=1e-12):
        self.tolerance = tolerance
        self.failed_tests = []

    def validate_parallel_result(self,
                                system_config: dict,
                                parallel_result: ManyBodyResult,
                                reference_key: str) -> bool:
        """Compare parallel result against sequential reference."""

        # Load reference result
        ref_result = self.load_reference(reference_key)

        # Comprehensive numerical comparison
        validation_results = {
            "total_energy": self._compare_floats(
                parallel_result.return_result,
                ref_result.return_result
            ),
            "nbody_energies": self._compare_nbody_dict(
                parallel_result.properties.return_energy_body_dict,
                ref_result.properties.return_energy_body_dict
            ),
            "component_properties": self._compare_component_properties(
                parallel_result.component_properties,
                ref_result.component_properties
            )
        }

        # Check for any failures
        all_passed = all(validation_results.values())
        if not all_passed:
            self._log_failure(system_config, validation_results)

        return all_passed

    def _compare_floats(self, a: float, b: float) -> bool:
        """Ultra-strict floating point comparison."""
        if not math.isclose(a, b, rel_tol=self.tolerance, abs_tol=self.tolerance):
            print(f"FLOAT MISMATCH: {a} vs {b}, diff={abs(a-b):.2e}")
            return False
        return True

    def _compare_arrays(self, a: np.ndarray, b: np.ndarray) -> bool:
        """Ultra-strict array comparison."""
        if not np.allclose(a, b, rtol=self.tolerance, atol=self.tolerance):
            diff = np.abs(a - b)
            max_diff = np.max(diff)
            print(f"ARRAY MISMATCH: max_diff={max_diff:.2e}")
            print(f"Mismatched elements: {np.sum(diff > self.tolerance)}")
            return False
        return True
```

### 3. Continuous Regression Testing

#### Development-Time Testing
Every code change must pass regression tests before merge:

```python
# Run after every commit during parallel development
class DevelopmentRegressionSuite:
    def test_dependency_graph_preserves_results(self):
        """Test that new dependency graph produces identical iteration order."""
        # Compare iterate_molecules() vs iterate_molecules_by_level()
        assert_identical_molecule_iteration()

    def test_parallel_multiprocessing_exact_match(self):
        """Test multiprocessing execution matches sequential."""
        for test_case in self.core_test_cases:
            seq_result = compute_sequential(test_case)
            par_result = compute_parallel(test_case, mode="multiprocessing")
            assert_results_identical(seq_result, par_result, tolerance=1e-12)

    def test_parallel_mpi_exact_match(self):
        """Test MPI execution matches sequential (when MPI available)."""
        if not has_mpi():
            pytest.skip("MPI not available")
        # Similar to multiprocessing test
```

#### Automated CI/CD Integration

```yaml
# GitHub Actions workflow addition
- name: Parallel Regression Tests
  run: |
    # Run complete regression suite
    pytest qcmanybody/tests/test_parallel_regression.py -v

    # Run parallel execution tests with different worker counts
    pytest qcmanybody/tests/test_parallel_execution.py::test_worker_count_1 -v
    pytest qcmanybody/tests/test_parallel_execution.py::test_worker_count_2 -v
    pytest qcmanybody/tests/test_parallel_execution.py::test_worker_count_4 -v

    # Validate numerical accuracy
    python scripts/validate_parallel_accuracy.py --tolerance 1e-12
```

### 4. Multi-Dimensional Validation Matrix

Test combinations of:

| Dimension | Values | Purpose |
|-----------|--------|---------|
| **System Size** | 2-8 fragments | Test scalability |
| **BSSE Treatment** | cp, nocp, vmfc | Test all correction types |
| **Max N-body** | 2, 3, 4 | Test dependency chains |
| **Driver** | energy, gradient, hessian | Test all property types |
| **Levels** | single, multi-level | Test complex level scenarios |
| **Workers** | 1, 2, 4, 8 | Test different parallelization |
| **QC Program** | Psi4, NWChem, CFOUR | Test all supported programs |

```python
@pytest.mark.parametrize("fragments", [2, 3, 4, 6])
@pytest.mark.parametrize("bsse_type", ["cp", "nocp", "vmfc"])
@pytest.mark.parametrize("max_nbody", [2, 3, 4])
@pytest.mark.parametrize("workers", [1, 2, 4])
def test_parallel_regression_matrix(fragments, bsse_type, max_nbody, workers):
    """Comprehensive regression test matrix."""

    # Skip invalid combinations
    if max_nbody > fragments:
        pytest.skip(f"max_nbody ({max_nbody}) > fragments ({fragments})")

    # Generate test system
    system = generate_test_system(fragments)
    config = {"max_nbody": max_nbody, "bsse_type": bsse_type}

    # Sequential reference
    seq_result = compute_sequential(system, config)

    # Parallel execution
    par_result = compute_parallel(system, config, workers=workers)

    # Strict validation
    assert_numerical_identity(seq_result, par_result, tolerance=1e-12)
```

### 5. Performance-Aware Regression Testing

Monitor that parallel improvements don't sacrifice correctness:

```python
class PerformanceRegressionTester:
    def test_parallel_correctness_vs_performance(self):
        """Ensure performance improvements maintain correctness."""

        test_cases = self.get_performance_test_cases()

        for case in test_cases:
            # Get reference result (high accuracy)
            ref_result = compute_with_tight_convergence(case)

            # Test parallel execution (normal accuracy)
            par_result = compute_parallel(case, workers=4)

            # Results should be identical to working precision
            assert_results_match(ref_result, par_result, tolerance=1e-10)

            # Performance should improve
            assert par_result.timing < ref_result.timing * 0.8  # At least 20% faster
```

### 6. Error Condition Testing

Ensure parallel execution handles errors consistently:

```python
class ErrorConditionTester:
    def test_parallel_error_consistency(self):
        """Test that parallel and sequential execution fail the same way."""

        error_cases = [
            {"invalid_geometry": "overlapping_atoms"},
            {"invalid_basis": "nonexistent_basis_set"},
            {"convergence_failure": "difficult_scf"},
        ]

        for error_case in error_cases:
            # Both should raise the same exception type
            with pytest.raises(Exception) as seq_exc:
                compute_sequential(error_case)

            with pytest.raises(Exception) as par_exc:
                compute_parallel(error_case, workers=2)

            # Exception types should match
            assert type(seq_exc.value) == type(par_exc.value)
```

## Implementation Recommendations

### Phase 1: Foundation Testing
1. **Generate Reference Dataset**: Create comprehensive golden references
2. **Basic Validation Framework**: Implement numerical comparison utilities
3. **Dependency Graph Testing**: Validate new iteration methods

### Phase 2: Parallel Execution Testing
1. **Worker Count Matrix**: Test 1, 2, 4, 8 workers systematically
2. **Load Balancing Validation**: Ensure different strategies give identical results
3. **Error Handling Tests**: Validate error conditions

### Phase 3: Integration Testing
1. **QC Program Matrix**: Test all supported QC programs
2. **Complex Scenarios**: Multi-level calculations, embedding charges
3. **Performance Regression**: Monitor performance while maintaining correctness

### Phase 4: Production Readiness
1. **Stress Testing**: Large systems, edge cases
2. **Long-term Stability**: Extended test runs
3. **User Acceptance**: Real-world scenario validation

## Test Data Management

```python
# Organize reference data systematically
reference_data/
├── he_systems/
│   ├── he2_cp_energy.json
│   ├── he4_vmfc_gradient.json
│   └── ...
├── water_clusters/
│   ├── h2o3_multilevel_energy.json
│   └── ...
└── validation_matrix/
    ├── system_2frag_cp_2body.json
    └── ...
```

## Success Criteria

1. **Zero Numerical Differences**: Parallel results must be identical to 1e-12 tolerance
2. **100% Test Coverage**: All existing tests must pass with parallel execution
3. **Performance Validation**: Parallel execution must be faster while maintaining accuracy
4. **Error Consistency**: Parallel and sequential execution must fail identically
5. **CI Integration**: All tests must pass in automated environment

This comprehensive regression testing strategy ensures that parallel development maintains the mathematical correctness that is absolutely critical for quantum chemistry calculations.

<function_calls>
<invoke name="TodoWrite">
<parameter name="todos">[{"content": "Analyze existing QCManyBody test suite structure", "status": "completed", "activeForm": "Analyzing existing QCManyBody test suite structure"}, {"content": "Design regression testing strategy for parallel development", "status": "completed", "activeForm": "Designing regression testing strategy for parallel development"}, {"content": "Create reference result validation framework", "status": "completed", "activeForm": "Creating reference result validation framework"}, {"content": "Design continuous integration testing approach", "status": "in_progress", "activeForm": "Designing continuous integration testing approach"}, {"content": "Document testing recommendations for the project", "status": "pending", "activeForm": "Documenting testing recommendations for the project"}]