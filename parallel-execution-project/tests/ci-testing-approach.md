# Continuous Integration Testing Approach

## Overview

This document outlines the CI/CD strategy for ensuring parallel execution development maintains correctness throughout the development process, with automated testing at every code change.

## Current QCManyBody CI Analysis

Based on the existing test infrastructure:
- Uses GitHub Actions (`.github/workflows/ci.yml`)
- Pytest with markers for external dependencies (`@pytest.mark.psi4`, `@pytest.mark.qcengine`)
- Test categories: unit tests, static-data regression, end-to-end with QC programs
- Supports multiple Python versions and platforms

## Enhanced CI Strategy for Parallel Development

### 1. Multi-Stage CI Pipeline

```yaml
# Enhanced GitHub Actions workflow
name: Parallel Execution CI

on: [push, pull_request]

jobs:
  # Stage 1: Fast validation (< 5 minutes)
  quick-validation:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, 3.10, 3.11, 3.12]
    steps:
      - uses: actions/checkout@v3
      - name: Set up Python
        uses: actions/setup-python@v3
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install dependencies
        run: |
          pip install -e .[tests]
          pip install zstandard  # For static regression tests

      - name: Unit tests
        run: pytest qcmanybody/tests/test_utils.py -v

      - name: Dependency graph tests
        run: pytest qcmanybody/tests/test_dependency.py -v

      - name: Static regression tests
        run: pytest qcmanybody/tests/test_core_*.py -v

  # Stage 2: Parallel execution validation (< 15 minutes)
  parallel-validation:
    needs: quick-validation
    runs-on: ubuntu-latest
    strategy:
      matrix:
        workers: [1, 2, 4]
    steps:
      - uses: actions/checkout@v3
      - name: Set up Python 3.10
        uses: actions/setup-python@v3
        with:
          python-version: "3.10"

      - name: Install test dependencies
        run: |
          pip install -e .[tests,parallel]
          pip install pytest-xdist  # For parallel test execution

      - name: Parallel regression tests
        run: |
          pytest qcmanybody/tests/test_parallel_regression.py \
            --workers=${{ matrix.workers }} -v

      - name: Numerical accuracy validation
        run: |
          python scripts/validate_parallel_accuracy.py \
            --tolerance=1e-12 --workers=${{ matrix.workers }}

  # Stage 3: QC program integration (< 30 minutes)
  qc-integration:
    needs: parallel-validation
    runs-on: ubuntu-latest
    strategy:
      matrix:
        qc-program: [psi4, nwchem]
        parallel-mode: [sequential, multiprocessing]
    steps:
      - uses: actions/checkout@v3
      - name: Set up conda environment
        uses: conda-incubator/setup-miniconda@v2
        with:
          environment-file: .github/envs/test-${{ matrix.qc-program }}.yml

      - name: Install QCManyBody
        run: pip install -e .[tests,parallel]
        shell: bash -l {0}

      - name: QC integration tests
        run: |
          pytest qcmanybody/tests/test_computer_*.py \
            -m "${{ matrix.qc-program }}" \
            --parallel-mode=${{ matrix.parallel-mode }} -v
        shell: bash -l {0}

      - name: Example tests
        run: |
          pytest qcmanybody/tests/test_examples.py \
            -m "${{ matrix.qc-program }}" \
            --parallel-mode=${{ matrix.parallel-mode }} -v
        shell: bash -l {0}

  # Stage 4: Performance regression (runs on main branch only)
  performance-regression:
    if: github.ref == 'refs/heads/main'
    needs: qc-integration
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
        with:
          fetch-depth: 0  # Need history for comparison

      - name: Set up conda environment
        uses: conda-incubator/setup-miniconda@v2
        with:
          environment-file: .github/envs/benchmark.yml

      - name: Performance benchmarks
        run: |
          python scripts/benchmark_parallel_execution.py \
            --compare-baseline --report-regression
        shell: bash -l {0}

      - name: Upload performance data
        uses: actions/upload-artifact@v3
        with:
          name: performance-data
          path: benchmarks/results/
```

### 2. Custom Pytest Markers and Configuration

```python
# pytest configuration for parallel testing
# pyproject.toml additions
[tool.pytest.ini_options]
markers = [
    # Existing markers
    "addon: tests require external non-required software",
    "qcengine: tests using QCEngine software; skip if unavailable",
    "psi4: tests using Psi4 software; skip if unavailable",
    "nwchem: tests using classic NWChem software; skip if unavailable",

    # New parallel execution markers
    "parallel: tests for parallel execution functionality",
    "regression: tests comparing parallel vs sequential results",
    "performance: tests measuring parallel execution performance",
    "slow: tests that take >30 seconds to run",
    "workers_1: tests with 1 worker (sequential mode)",
    "workers_2: tests with 2 workers",
    "workers_4: tests with 4 workers",
    "workers_8: tests with 8 workers",
    "mpi: tests requiring MPI installation",
    "memory_intensive: tests requiring >4GB RAM",
]

# Pytest fixtures for parallel testing
@pytest.fixture(scope="session")
def reference_data_cache():
    """Cache reference results for regression testing."""
    return ReferenceDataCache()

@pytest.fixture(params=[1, 2, 4])
def worker_count(request):
    """Parameterize tests across different worker counts."""
    return request.param

@pytest.fixture
def parallel_config(worker_count):
    """Create parallel configuration for testing."""
    return ParallelConfig(
        mode=ParallelMode.MULTIPROCESSING,
        max_workers=worker_count,
        timeout_per_fragment=300
    )
```

### 3. Automated Reference Data Management

```python
# CI script for managing reference data
class CIReferenceManager:
    def __init__(self):
        self.reference_dir = "qcmanybody/tests/reference_data"
        self.current_version = self.get_git_commit()

    def validate_references_current(self):
        """Ensure reference data is up to date."""
        ref_version = self.load_reference_version()

        if ref_version != self.current_version:
            # Check if reference update is needed
            if self.sequential_code_changed():
                raise ValueError(
                    "Sequential code changed but reference data not updated. "
                    "Run: python scripts/update_reference_data.py"
                )

    def generate_references_if_missing(self):
        """Auto-generate reference data for new test cases."""
        for test_case in self.discover_test_cases():
            ref_file = f"{self.reference_dir}/{test_case.name}.json"
            if not os.path.exists(ref_file):
                print(f"Generating reference data for {test_case.name}")
                reference_result = self.compute_reference(test_case)
                self.save_reference(ref_file, reference_result)
```

### 4. Parallel-Specific Test Categories

#### A. Regression Tests (Always Run)
```python
@pytest.mark.regression
class TestParallelRegression:
    """Tests ensuring parallel execution matches sequential."""

    @pytest.mark.parametrize("system", ["he2", "h2o3", "ar4"])
    @pytest.mark.parametrize("bsse", ["cp", "nocp", "vmfc"])
    def test_parallel_vs_sequential_exact_match(self, system, bsse, worker_count):
        """Core regression test - must pass always."""
        seq_result = compute_sequential(system, bsse_type=bsse)
        par_result = compute_parallel(system, bsse_type=bsse, workers=worker_count)

        assert_numerical_identity(seq_result, par_result, tolerance=1e-12)

@pytest.mark.regression
@pytest.mark.slow
class TestComprehensiveRegression:
    """Comprehensive regression tests for complex scenarios."""

    def test_multilevel_parallel_regression(self):
        """Test multi-level calculations in parallel."""
        # Complex multi-level test cases

    def test_embedding_charges_parallel_regression(self):
        """Test embedding charges in parallel execution."""
        # Embedding charge scenarios
```

#### B. Performance Tests (Selective Running)
```python
@pytest.mark.performance
class TestParallelPerformance:
    """Performance validation tests."""

    @pytest.mark.slow
    def test_parallel_speedup_validation(self, benchmark):
        """Validate parallel execution provides speedup."""
        system = create_large_test_system(fragments=6)

        # Benchmark sequential
        seq_time = benchmark.pedantic(
            compute_sequential, args=(system,), rounds=3
        )

        # Benchmark parallel
        par_time = benchmark.pedantic(
            compute_parallel, args=(system,), kwargs={"workers": 4}, rounds=3
        )

        speedup = seq_time / par_time
        assert speedup > 1.5  # At least 1.5x speedup
```

#### C. Infrastructure Tests (Development Only)
```python
@pytest.mark.parallel
class TestParallelInfrastructure:
    """Test parallel execution infrastructure."""

    def test_worker_process_management(self):
        """Test worker processes start and stop correctly."""

    def test_memory_limit_enforcement(self):
        """Test memory limits are enforced per worker."""

    def test_error_handling_parallel(self):
        """Test error handling in parallel execution."""
```

### 5. Environment-Specific CI Configurations

#### Development Environment
```yaml
# .github/envs/test-dev.yml
name: qcmanybody-dev-test
channels:
  - conda-forge
dependencies:
  - python=3.10
  - pytest
  - pytest-xdist
  - pytest-benchmark
  - zstandard
  - numpy
  - pydantic>=1.10,<3
  - qcelemental
  - pip:
    - -e .[parallel]
```

#### QC Program Environments
```yaml
# .github/envs/test-psi4.yml
name: qcmanybody-psi4-test
channels:
  - conda-forge
  - psi4
dependencies:
  - python=3.10
  - psi4=1.9.1
  - qcengine
  - pytest
  - pip:
    - -e .[parallel]

# .github/envs/test-nwchem.yml
name: qcmanybody-nwchem-test
channels:
  - conda-forge
dependencies:
  - python=3.10
  - nwchem
  - qcengine
  - pytest
  - pip:
    - -e .[parallel]
```

### 6. CI Failure Handling and Reporting

```python
# CI failure analysis script
class CIFailureAnalyzer:
    def analyze_test_failure(self, test_output):
        """Analyze test failures and provide debugging info."""

        if "numerical mismatch" in test_output.lower():
            return self.analyze_numerical_failure(test_output)
        elif "worker process" in test_output.lower():
            return self.analyze_worker_failure(test_output)
        elif "timeout" in test_output.lower():
            return self.analyze_timeout_failure(test_output)

    def generate_failure_report(self, failure_info):
        """Generate detailed failure report for developers."""
        report = {
            "failure_type": failure_info.type,
            "system_info": self.get_system_info(),
            "reproduction_steps": failure_info.reproduction_steps,
            "debugging_hints": failure_info.debugging_hints
        }
        return report
```

### 7. Performance Regression Detection

```python
# Performance regression detection
class PerformanceMonitor:
    def __init__(self):
        self.baseline_file = "benchmarks/baseline_performance.json"
        self.regression_threshold = 0.15  # 15% slowdown threshold

    def check_performance_regression(self, current_results):
        """Check if current performance has regressed."""
        baseline = self.load_baseline()

        for test_case, current_time in current_results.items():
            baseline_time = baseline.get(test_case)

            if baseline_time and current_time > baseline_time * (1 + self.regression_threshold):
                raise PerformanceRegressionError(
                    f"Performance regression in {test_case}: "
                    f"{current_time:.2f}s vs baseline {baseline_time:.2f}s"
                )

    def update_baseline(self, new_results):
        """Update baseline performance data."""
        # Only update if performance improved
        baseline = self.load_baseline()
        updated = False

        for test_case, new_time in new_results.items():
            if test_case not in baseline or new_time < baseline[test_case]:
                baseline[test_case] = new_time
                updated = True

        if updated:
            self.save_baseline(baseline)
```

### 8. CI Success Criteria

#### Every Commit Must Pass:
- [ ] All unit tests pass
- [ ] All static regression tests pass
- [ ] New dependency graph tests pass
- [ ] Basic parallel execution tests pass (1, 2 workers)

#### Pull Request Must Pass:
- [ ] All commit-level tests pass
- [ ] Parallel regression tests pass (1, 2, 4 workers)
- [ ] At least one QC program integration test passes
- [ ] No performance regression detected

#### Main Branch Must Pass:
- [ ] All pull request tests pass
- [ ] All QC program integration tests pass
- [ ] Full performance benchmark suite passes
- [ ] Long-running stability tests pass

This comprehensive CI approach ensures that parallel development maintains numerical correctness while providing rapid feedback to developers and preventing regression bugs from reaching the main branch.

<function_calls>
<invoke name="TodoWrite">
<parameter name="todos">[{"content": "Analyze existing QCManyBody test suite structure", "status": "completed", "activeForm": "Analyzing existing QCManyBody test suite structure"}, {"content": "Design regression testing strategy for parallel development", "status": "completed", "activeForm": "Designing regression testing strategy for parallel development"}, {"content": "Create reference result validation framework", "status": "completed", "activeForm": "Creating reference result validation framework"}, {"content": "Design continuous integration testing approach", "status": "completed", "activeForm": "Designing continuous integration testing approach"}, {"content": "Document testing recommendations for the project", "status": "in_progress", "activeForm": "Documenting testing recommendations for the project"}]