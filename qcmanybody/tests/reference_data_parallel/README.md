# QCManyBody Parallel Execution Reference Dataset

This directory contains the **Golden Reference Dataset** for QCManyBody parallel execution validation. This dataset represents the absolute truth for regression testing during parallel development.

## 🎯 Purpose

The reference dataset ensures that parallel execution produces **mathematically identical** results to sequential execution with ultra-strict precision (1e-12 tolerance).

## 📁 Dataset Structure

```
reference_data_parallel/
├── reference_data_v1.0.json.zst    # Compressed main dataset
├── reference_metadata.json         # Dataset metadata and generation info
├── reference_index.json           # Test case index and navigation
└── README.md                      # This file
```

## 🧪 Test Case Coverage

The dataset systematically covers:

- **Systems**: ne2 (simple dimer), h2o3 (water trimer)
- **BSSE Types**: cp, nocp, vmfc
- **N-body Levels**: 2-body, 3-body
- **Drivers**: energy, gradient
- **Level Combinations**: single-level (SCF), single-level (MP2), mixed (SCF→MP2)

**Total Expected Cases**: ~72 systematic combinations

## 🔬 Numerical Precision

- **Generation Precision**: 1e-14+ (working precision)
- **Validation Tolerance**: 1e-12 (ultra-strict)
- **Sequential Code**: Unmodified QCManyBody main branch

## 🚀 Usage

### Loading Reference Data
```python
from qcmanybody.testing import ReferenceDataLoader

loader = ReferenceDataLoader()
reference_case = loader.get_reference_case("ne2_cp_energy_maxnb2_single_scf")
```

### Validating Results
```python
from qcmanybody.testing import ParallelRegressionTester

tester = ParallelRegressionTester(tolerance=1e-12)
report = tester.validate_result(parallel_result, "test_case_id")
print(report.generate_summary())
```

### Batch Validation
```python
batch_report = tester.validate_batch(all_parallel_results, reference_keys)
print(batch_report.generate_summary())
```

## 📊 Dataset Generation

Generate the dataset:
```bash
python scripts/generate_reference_data.py
```

Validate dataset integrity:
```bash
python scripts/validate_reference_data.py --all
```

## ⚠️ Critical Requirements

1. **Never modify reference data** once generated and validated
2. **Use only unmodified sequential code** for generation
3. **Maintain 1e-12 tolerance** for all validations
4. **Version control** all reference data changes

## 🔍 Data Format

Each reference case contains:
```json
{
  "test_case_id": "system_bsse_driver_maxnbX_levels",
  "system": {
    "name": "system_name",
    "description": "Human-readable description",
    "molecule": {...}
  },
  "configuration": {
    "bsse_type": "cp|nocp|vmfc",
    "max_nbody": 2|3|4,
    "driver": "energy|gradient|hessian",
    "levels": {1: "method", 2: "method", ...}
  },
  "component_results": {...},
  "final_results": {
    "ret_energy": -123.456789,
    "energy_body_dict": {...},
    "results": {...}
  },
  "numerical_precision": 1e-14
}
```

## 🧪 Validation Framework

The validation framework provides:
- **Ultra-strict numerical comparison** (1e-12 tolerance)
- **Detailed failure reporting** with exact difference locations
- **Element-wise array analysis** for gradients/hessians
- **Performance metrics** and timing information
- **Batch validation** for comprehensive testing

## 📈 Quality Assurance

- **Automatic integrity checks** verify dataset completeness
- **Precision validation** ensures numerical requirements
- **Cross-validation** with existing test suite
- **Version tracking** maintains dataset provenance

## 🔧 Troubleshooting

### Missing Reference Data
```bash
# Generate the dataset
python scripts/generate_reference_data.py

# Verify generation
python scripts/validate_reference_data.py --verify-integrity
```

### Validation Failures
```bash
# Check tolerance settings
python -c "from qcmanybody.testing import ParallelRegressionTester; print(ParallelRegressionTester().tolerance)"

# Generate detailed report
python scripts/validate_reference_data.py --generate-report --output validation_report.json
```

### Dataset Corruption
```bash
# Verify checksums and structure
python scripts/validate_reference_data.py --all

# Regenerate if necessary (CAUTION: affects all dependent testing)
python scripts/generate_reference_data.py --force-regenerate
```

## 📚 Related Documentation

- **Phase 0 Task**: `parallel-execution-project/tasks/phase-0/T0-001-reference-dataset.md`
- **Validation Framework**: `parallel-execution-project/tasks/phase-0/T0-002-validation-framework.md`
- **Testing Strategy**: `parallel-execution-project/tests/regression-testing-strategy.md`
- **API Documentation**: `qcmanybody/testing/__init__.py`

---

**⚡ Remember**: This dataset is the foundation for all parallel execution validation. Any changes must be carefully considered and thoroughly tested!