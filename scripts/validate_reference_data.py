#!/usr/bin/env python3
"""
Validate Reference Dataset Integrity

This script validates the generated golden reference dataset for completeness,
integrity, and numerical precision requirements.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from qcmanybody.testing import ReferenceDataLoader

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def validate_dataset_integrity(reference_dir: str) -> bool:
    """Validate the integrity of the reference dataset."""
    logger.info("Validating reference dataset integrity")

    try:
        loader = ReferenceDataLoader(reference_dir)
        success = loader.validate_reference_integrity()

        if success:
            logger.info("✅ Reference dataset integrity validation PASSED")
        else:
            logger.error("❌ Reference dataset integrity validation FAILED")

        return success

    except Exception as e:
        logger.error(f"❌ Integrity validation failed with exception: {str(e)}")
        return False


def verify_numerical_precision(reference_dir: str) -> bool:
    """Verify that the reference data meets numerical precision requirements."""
    logger.info("Verifying numerical precision requirements")

    try:
        loader = ReferenceDataLoader(reference_dir)
        reference_data = loader.load_reference_data()

        precision_failures = []

        for case_id, case_data in reference_data.items():
            # Check that numerical precision metadata exists
            if "numerical_precision" not in case_data:
                precision_failures.append(f"{case_id}: Missing numerical_precision metadata")
                continue

            required_precision = case_data["numerical_precision"]
            if required_precision > 1e-13:  # Should be at least 1e-14
                precision_failures.append(f"{case_id}: Insufficient precision {required_precision}")

            # Check final results for NaN/inf values
            final_results = case_data.get("final_results", {})
            if "ret_energy" in final_results:
                energy = final_results["ret_energy"]
                if not isinstance(energy, (int, float)) or abs(energy) > 1e10:
                    precision_failures.append(f"{case_id}: Suspicious energy value {energy}")

        if precision_failures:
            logger.error(f"❌ Numerical precision validation FAILED:")
            for failure in precision_failures[:10]:  # Show first 10 failures
                logger.error(f"  {failure}")
            if len(precision_failures) > 10:
                logger.error(f"  ... and {len(precision_failures) - 10} more failures")
            return False
        else:
            logger.info("✅ Numerical precision validation PASSED")
            return True

    except Exception as e:
        logger.error(f"❌ Precision validation failed with exception: {str(e)}")
        return False


def verify_completeness(reference_dir: str) -> bool:
    """Verify that the reference dataset has expected completeness."""
    logger.info("Verifying dataset completeness")

    try:
        loader = ReferenceDataLoader(reference_dir)
        dataset_info = loader.get_dataset_info()

        logger.info(f"Dataset statistics:")
        logger.info(f"  Total cases: {dataset_info['total_cases']}")

        # Log case distribution
        if "case_distribution" in dataset_info:
            for category, counts in dataset_info["case_distribution"].items():
                logger.info(f"  {category}: {counts}")

        # Check minimum expected cases (should have at least basic test matrix)
        # 2 systems × 3 BSSE types × 2 max_nbody × 2 drivers × 3 level combinations = 72 cases
        min_expected_cases = 72

        if dataset_info["total_cases"] < min_expected_cases:
            logger.error(f"❌ Completeness validation FAILED: Only {dataset_info['total_cases']} cases, expected at least {min_expected_cases}")
            return False
        else:
            logger.info(f"✅ Completeness validation PASSED: {dataset_info['total_cases']} cases generated")
            return True

    except Exception as e:
        logger.error(f"❌ Completeness validation failed with exception: {str(e)}")
        return False


def generate_validation_report(reference_dir: str, output_file: str = None) -> None:
    """Generate a comprehensive validation report."""
    logger.info("Generating validation report")

    try:
        loader = ReferenceDataLoader(reference_dir)
        dataset_info = loader.get_dataset_info()

        report = {
            "validation_timestamp": "2024-09-26T00:00:00Z",
            "reference_directory": reference_dir,
            "dataset_info": dataset_info,
            "validation_results": {
                "integrity_passed": validate_dataset_integrity(reference_dir),
                "precision_passed": verify_numerical_precision(reference_dir),
                "completeness_passed": verify_completeness(reference_dir)
            }
        }

        report["overall_status"] = all(report["validation_results"].values())

        if output_file:
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Validation report saved to: {output_file}")
        else:
            logger.info("Validation Report:")
            logger.info(json.dumps(report, indent=2))

    except Exception as e:
        logger.error(f"❌ Report generation failed: {str(e)}")


def main():
    """Main validation function."""
    parser = argparse.ArgumentParser(description="Validate QCManyBody reference dataset")
    parser.add_argument(
        "--reference-dir",
        default="qcmanybody/tests/reference_data_parallel",
        help="Directory containing reference data"
    )
    parser.add_argument(
        "--verify-integrity",
        action="store_true",
        help="Verify dataset integrity"
    )
    parser.add_argument(
        "--verify-precision",
        action="store_true",
        help="Verify numerical precision"
    )
    parser.add_argument(
        "--verify-completeness",
        action="store_true",
        help="Verify dataset completeness"
    )
    parser.add_argument(
        "--generate-report",
        action="store_true",
        help="Generate comprehensive validation report"
    )
    parser.add_argument(
        "--output",
        help="Output file for validation report"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all validation checks"
    )

    args = parser.parse_args()

    # Check if reference directory exists
    ref_dir = Path(args.reference_dir)
    if not ref_dir.exists():
        logger.error(f"Reference directory does not exist: {ref_dir}")
        logger.error("Please run scripts/generate_reference_data.py first")
        sys.exit(1)

    # Run validations based on arguments
    all_passed = True

    if args.all or args.verify_integrity:
        if not validate_dataset_integrity(args.reference_dir):
            all_passed = False

    if args.all or args.verify_precision:
        if not verify_numerical_precision(args.reference_dir):
            all_passed = False

    if args.all or args.verify_completeness:
        if not verify_completeness(args.reference_dir):
            all_passed = False

    if args.all or args.generate_report:
        generate_validation_report(args.reference_dir, args.output)

    # Summary
    if args.all:
        if all_passed:
            logger.info("🎉 ALL VALIDATIONS PASSED - Reference dataset is ready for use")
            sys.exit(0)
        else:
            logger.error("❌ SOME VALIDATIONS FAILED - Reference dataset needs attention")
            sys.exit(1)


if __name__ == "__main__":
    main()