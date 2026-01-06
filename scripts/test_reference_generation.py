#!/usr/bin/env python3
"""
Test Reference Generation System

This script tests the reference generation system with a minimal test case
to ensure the infrastructure works before running the full dataset generation.
"""

import logging
import sys
import time
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from qcmanybody.testing import ReferenceDataLoader, ParallelRegressionTester

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def test_validation_framework():
    """Test the validation framework with known values."""
    logger.info("Testing validation framework")

    tester = ParallelRegressionTester(tolerance=1e-12)

    # Test float comparison
    result = tester.compare_floats(1.0, 1.0)
    assert result.passed, "Identical floats should pass"

    result = tester.compare_floats(1.0, 1.0 + 1e-13)
    assert result.passed, "Values within tolerance should pass"

    result = tester.compare_floats(1.0, 1.0 + 1e-11)
    assert not result.passed, "Values outside tolerance should fail"

    # Test array comparison
    import numpy as np

    a = np.array([1.0, 2.0, 3.0])
    b = np.array([1.0, 2.0, 3.0])
    result = tester.compare_arrays(a, b)
    assert result.passed, "Identical arrays should pass"

    b = np.array([1.0, 2.0, 3.0 + 1e-11])
    result = tester.compare_arrays(a, b)
    assert not result.passed, "Arrays outside tolerance should fail"

    logger.info("✅ Validation framework tests passed")


def test_existing_validation():
    """Test against existing QCManyBody validation to ensure compatibility."""
    logger.info("Testing compatibility with existing validation")

    # Import existing validation function
    from qcmanybody.tests.utils import compare

    # Test that our stricter validation catches differences the old one misses
    try:
        # This should pass with old validation (1e-7 tolerance)
        compare(1.0, 1.0 + 5e-8)
        logger.info("Old validation passes 5e-8 difference as expected")
    except:
        logger.error("Old validation unexpectedly failed")

    # Our validation should fail
    tester = ParallelRegressionTester(tolerance=1e-12)
    result = tester.compare_floats(1.0, 1.0 + 5e-8)
    assert not result.passed, "New validation should be stricter"

    logger.info("✅ Validation compatibility tests passed")


def test_data_loading_structure():
    """Test the data loading structure without requiring actual data."""
    logger.info("Testing data loading structure")

    # Test that loader can be instantiated
    try:
        loader = ReferenceDataLoader("nonexistent_directory")
        logger.info("✅ ReferenceDataLoader instantiated successfully")
    except Exception as e:
        logger.error(f"❌ Failed to instantiate ReferenceDataLoader: {e}")
        return False

    # Test error handling for missing data
    try:
        loader.load_reference_data()
        logger.error("❌ Should have failed on missing data")
        return False
    except FileNotFoundError:
        logger.info("✅ Proper error handling for missing data")

    return True


def check_dependencies():
    """Check that required dependencies are available."""
    logger.info("Checking dependencies")

    required_modules = [
        "qcelemental",
        "qcmanybody",
        "numpy",
        "zstandard"
    ]

    missing_modules = []
    for module in required_modules:
        try:
            __import__(module)
            logger.info(f"✅ {module} available")
        except ImportError:
            missing_modules.append(module)
            logger.error(f"❌ {module} missing")

    if missing_modules:
        logger.error(f"Missing required modules: {missing_modules}")
        return False

    # Check optional QC modules (needed for actual reference generation)
    try:
        import qcengine
        logger.info("✅ qcengine available")
    except ImportError:
        logger.warning("⚠️  qcengine not available - needed for reference generation")

    try:
        import psi4
        logger.info("✅ psi4 available")
    except ImportError:
        logger.warning("⚠️  psi4 not available - needed for QC calculations")

    return True


def main():
    """Run all tests."""
    logger.info("Starting reference generation system tests")

    all_tests_passed = True

    # Check dependencies
    if not check_dependencies():
        logger.error("❌ Dependency check failed")
        all_tests_passed = False

    # Test validation framework
    try:
        test_validation_framework()
    except Exception as e:
        logger.error(f"❌ Validation framework test failed: {e}")
        all_tests_passed = False

    # Test compatibility
    try:
        test_existing_validation()
    except Exception as e:
        logger.error(f"❌ Compatibility test failed: {e}")
        all_tests_passed = False

    # Test data loading
    try:
        test_data_loading_structure()
    except Exception as e:
        logger.error(f"❌ Data loading test failed: {e}")
        all_tests_passed = False

    # Summary
    if all_tests_passed:
        logger.info("🎉 ALL TESTS PASSED - Reference generation system is ready")
        logger.info("You can now run: python scripts/generate_reference_data.py")
        return True
    else:
        logger.error("❌ SOME TESTS FAILED - Fix issues before proceeding")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)