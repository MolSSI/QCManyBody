#!/usr/bin/env python3
"""Debug validation issue."""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from qcmanybody.testing import ReferenceDataLoader, ParallelRegressionTester

def debug_validation():
    """Debug the validation failure."""

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Create the same mock data
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

        # Save data
        import zstandard
        reference_file = temp_path / "reference_data_v1.0.json.zst"
        with zstandard.open(reference_file, "wt") as f:
            json.dump(mock_reference_data, f, indent=2)

        # Test step by step
        loader = ReferenceDataLoader(temp_path)
        case_data = loader.get_reference_case("test_case_minimal")

        print("Reference case data:")
        print(json.dumps(case_data, indent=2))

        tester = ParallelRegressionTester(tolerance=1e-12, reference_data_path=temp_path)

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

        print("\nParallel result:")
        print(json.dumps(mock_parallel_result, indent=2))

        report = tester.validate_result(mock_parallel_result, "test_case_minimal")

        print(f"\nValidation result: {report.passed}")
        print(f"Error count: {len(report.error_details)}")

        for i, error in enumerate(report.error_details):
            print(f"\nError {i+1}:")
            print(f"  Location: {error.location}")
            print(f"  Expected: {error.expected_value}")
            print(f"  Actual: {error.actual_value}")
            print(f"  Type: {error.error_type}")
            print(f"  Description: {error.description}")

if __name__ == "__main__":
    debug_validation()