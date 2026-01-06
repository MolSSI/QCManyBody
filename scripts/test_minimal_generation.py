#!/usr/bin/env python3
"""
Test minimal reference generation to verify end-to-end functionality.
"""

import json
import logging
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from qcmanybody.testing import ReferenceDataLoader, ParallelRegressionTester

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def test_minimal_reference_generation():
    """Test reference generation with a mock dataset."""
    logger.info("Testing minimal reference generation")

    # Create a temporary directory for testing
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Create mock reference data structure
        mock_reference_data = {
            "test_case_minimal": {
                "test_case_id": "test_case_minimal",
                "system": {
                    "name": "mock_system",
                    "description": "Mock system for testing",
                    "molecule": {"symbols": ["H", "H"], "geometry": [0, 0, 0, 0, 0, 1]}
                },
                "configuration": {
                    "bsse_type": "nocp",
                    "max_nbody": 2,
                    "driver": "energy",
                    "levels": {1: "hf"},
                    "levels_name": "single_hf"
                },
                "final_results": {
                    "ret_energy": -1.123456789012345,
                    "energy_body_dict": {
                        "nocp": {1: -0.5, 2: -0.1}
                    },
                    "results": {
                        "NOCP-CORRECTED TOTAL ENERGY": -1.123456789012345,
                        "NOCP-CORRECTED INTERACTION ENERGY": -0.1
                    }
                },
                "numerical_precision": 1e-14,
                "timestamp": "2024-09-26T00:00:00Z"
            }
        }

        mock_metadata = {
            "generator_version": "1.0.0",
            "total_test_cases": 1,
            "generation_timestamp": "2024-09-26T00:00:00Z",
            "numerical_precision": 1e-14
        }

        mock_index = {
            "dataset_version": "1.0.0",
            "total_cases": 1,
            "test_case_index": {
                "test_case_minimal": {
                    "system": "mock_system",
                    "bsse_type": "nocp",
                    "max_nbody": 2,
                    "driver": "energy",
                    "levels": "single_hf"
                }
            }
        }

        # Save mock data using the same format as real system
        import zstandard

        reference_file = temp_path / "reference_data_v1.0.json.zst"
        with zstandard.open(reference_file, "wt") as f:
            json.dump(mock_reference_data, f, indent=2)

        metadata_file = temp_path / "reference_metadata.json"
        with open(metadata_file, "w") as f:
            json.dump(mock_metadata, f, indent=2)

        index_file = temp_path / "reference_index.json"
        with open(index_file, "w") as f:
            json.dump(mock_index, f, indent=2)

        logger.info(f"Created mock reference data in {temp_path}")

        # Test loading the data
        loader = ReferenceDataLoader(temp_path)

        # Test integrity validation
        try:
            integrity_ok = loader.validate_reference_integrity()
            assert integrity_ok, "Integrity validation should pass"
            logger.info("✅ Reference data integrity validation passed")
        except Exception as e:
            logger.error(f"❌ Integrity validation failed: {e}")
            return False

        # Test loading specific case
        try:
            case_data = loader.get_reference_case("test_case_minimal")
            assert case_data["test_case_id"] == "test_case_minimal"
            logger.info("✅ Reference case loading passed")
        except Exception as e:
            logger.error(f"❌ Case loading failed: {e}")
            return False

        # Test validation framework with this data
        try:
            tester = ParallelRegressionTester(tolerance=1e-12, reference_data_path=temp_path)

            # Create mock parallel result that should match
            mock_parallel_result = {
                "ret_energy": -1.123456789012345,
                "energy_body_dict": {
                    "nocp": {1: -0.5, 2: -0.1}
                },
                "results": {
                    "NOCP-CORRECTED TOTAL ENERGY": -1.123456789012345,
                    "NOCP-CORRECTED INTERACTION ENERGY": -0.1
                }
            }

            report = tester.validate_result(mock_parallel_result, "test_case_minimal")
            assert report.passed, f"Validation should pass, but got: {report.generate_summary()}"
            logger.info("✅ Validation framework with mock data passed")

            # Test with slightly different result (should fail)
            mock_parallel_result_bad = mock_parallel_result.copy()
            mock_parallel_result_bad["ret_energy"] = -1.123456789012345 + 1e-10  # Outside tolerance

            report_bad = tester.validate_result(mock_parallel_result_bad, "test_case_minimal")
            assert not report_bad.passed, "Validation should fail for out-of-tolerance result"
            logger.info("✅ Validation correctly detects differences outside tolerance")

        except Exception as e:
            logger.error(f"❌ Validation framework test failed: {e}")
            return False

        logger.info("✅ All end-to-end tests passed")
        return True


def test_validation_edge_cases():
    """Test validation framework edge cases."""
    logger.info("Testing validation edge cases")

    tester = ParallelRegressionTester(tolerance=1e-12)

    # Test handling of missing properties
    test_result = {"ret_energy": 1.0}
    reference_result = {"ret_energy": 1.0, "ret_gradient": [1.0, 2.0, 3.0]}

    # This should create a report with errors for missing gradient
    from qcmanybody.testing.validation_report import ValidationReport
    report = ValidationReport(
        test_id="test_missing_props",
        reference_key="test_missing_props",
        passed=True,
        tolerance_used=1e-12
    )

    try:
        tester._validate_final_results(test_result, reference_result, report)
        assert not report.passed, "Should fail when missing properties"
        assert len(report.error_details) > 0, "Should have error details"
        logger.info("✅ Correctly handles missing properties")
    except Exception as e:
        logger.error(f"❌ Missing property test failed: {e}")
        return False

    # Test array shape mismatch
    import numpy as np
    result = tester.compare_arrays(np.array([1, 2]), np.array([[1, 2]]))
    assert not result.passed, "Should fail on shape mismatch"
    logger.info("✅ Correctly detects array shape mismatches")

    return True


def main():
    """Run all end-to-end tests."""
    logger.info("Starting comprehensive end-to-end validation")

    tests_passed = 0
    total_tests = 2

    # Test 1: Minimal reference generation
    try:
        if test_minimal_reference_generation():
            tests_passed += 1
            logger.info("✅ Test 1 PASSED: Minimal reference generation")
        else:
            logger.error("❌ Test 1 FAILED: Minimal reference generation")
    except Exception as e:
        logger.error(f"❌ Test 1 FAILED with exception: {e}")

    # Test 2: Validation edge cases
    try:
        if test_validation_edge_cases():
            tests_passed += 1
            logger.info("✅ Test 2 PASSED: Validation edge cases")
        else:
            logger.error("❌ Test 2 FAILED: Validation edge cases")
    except Exception as e:
        logger.error(f"❌ Test 2 FAILED with exception: {e}")

    # Summary
    logger.info(f"End-to-end tests complete: {tests_passed}/{total_tests} passed")

    if tests_passed == total_tests:
        logger.info("🎉 ALL END-TO-END TESTS PASSED - System is fully functional")
        return True
    else:
        logger.error("❌ SOME END-TO-END TESTS FAILED - System needs fixes")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)