"""
Reference data loading utilities for accessing golden reference dataset.
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import zstandard


class ReferenceDataLoader:
    """Load and manage access to golden reference dataset."""

    def __init__(self, reference_dir: Union[str, Path] = None):
        """Initialize reference data loader.

        Args:
            reference_dir: Directory containing reference data files.
                         Defaults to qcmanybody/tests/reference_data_parallel/
        """
        if reference_dir is None:
            # Default to the standard reference data directory
            current_dir = Path(__file__).parent.parent
            reference_dir = current_dir / "tests" / "reference_data_parallel"

        self.reference_dir = Path(reference_dir)
        self._reference_data = None
        self._metadata = None
        self._index = None

    def load_reference_data(self, force_reload: bool = False) -> Dict[str, Any]:
        """Load the main reference dataset."""
        if self._reference_data is None or force_reload:
            reference_file = self.reference_dir / "reference_data_v1.0.json.zst"

            if not reference_file.exists():
                raise FileNotFoundError(
                    f"Reference data file not found: {reference_file}\n"
                    f"Please run scripts/generate_reference_data.py first"
                )

            with zstandard.open(reference_file, "rt") as f:
                self._reference_data = json.load(f)

        return self._reference_data

    def load_metadata(self, force_reload: bool = False) -> Dict[str, Any]:
        """Load reference dataset metadata."""
        if self._metadata is None or force_reload:
            metadata_file = self.reference_dir / "reference_metadata.json"

            if not metadata_file.exists():
                raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

            with open(metadata_file, "r") as f:
                self._metadata = json.load(f)

        return self._metadata

    def load_index(self, force_reload: bool = False) -> Dict[str, Any]:
        """Load reference dataset index."""
        if self._index is None or force_reload:
            index_file = self.reference_dir / "reference_index.json"

            if not index_file.exists():
                raise FileNotFoundError(f"Index file not found: {index_file}")

            with open(index_file, "r") as f:
                self._index = json.load(f)

        return self._index

    def get_reference_case(self, test_case_id: str) -> Dict[str, Any]:
        """Get reference data for a specific test case."""
        reference_data = self.load_reference_data()

        if test_case_id not in reference_data:
            available_cases = list(reference_data.keys())
            raise KeyError(
                f"Test case '{test_case_id}' not found in reference data.\n"
                f"Available cases: {available_cases[:5]}... (showing first 5)"
            )

        return reference_data[test_case_id]

    def list_test_cases(self, filter_by: Optional[Dict[str, Any]] = None) -> List[str]:
        """List available test cases, optionally filtered by criteria."""
        index = self.load_index()
        test_cases = list(index["test_case_index"].keys())

        if filter_by is None:
            return test_cases

        # Filter test cases based on criteria
        filtered_cases = []
        for case_id in test_cases:
            case_info = index["test_case_index"][case_id]
            match = True

            for key, value in filter_by.items():
                if key in case_info and case_info[key] != value:
                    match = False
                    break

            if match:
                filtered_cases.append(case_id)

        return filtered_cases

    def get_cases_by_system(self, system_name: str) -> List[str]:
        """Get all test cases for a specific system."""
        return self.list_test_cases(filter_by={"system": system_name})

    def get_cases_by_bsse_type(self, bsse_type: str) -> List[str]:
        """Get all test cases for a specific BSSE type."""
        return self.list_test_cases(filter_by={"bsse_type": bsse_type})

    def get_cases_by_driver(self, driver: str) -> List[str]:
        """Get all test cases for a specific driver (energy/gradient/hessian)."""
        return self.list_test_cases(filter_by={"driver": driver})

    def validate_reference_integrity(self) -> bool:
        """Validate the integrity of the reference dataset."""
        try:
            # Load all components
            reference_data = self.load_reference_data()
            metadata = self.load_metadata()
            index = self.load_index()

            # Check metadata consistency
            if len(reference_data) != metadata.get("total_test_cases", 0):
                raise ValueError("Metadata total_test_cases doesn't match actual data")

            if len(reference_data) != index.get("total_cases", 0):
                raise ValueError("Index total_cases doesn't match actual data")

            # Check that all index entries exist in reference data
            index_cases = set(index["test_case_index"].keys())
            reference_cases = set(reference_data.keys())

            if index_cases != reference_cases:
                missing_in_ref = index_cases - reference_cases
                missing_in_index = reference_cases - index_cases
                raise ValueError(
                    f"Index and reference data mismatch.\n"
                    f"Missing in reference: {missing_in_ref}\n"
                    f"Missing in index: {missing_in_index}"
                )

            # Validate each reference case has required structure
            required_keys = ["test_case_id", "system", "configuration", "final_results"]
            for case_id, case_data in reference_data.items():
                for key in required_keys:
                    if key not in case_data:
                        raise ValueError(f"Case {case_id} missing required key: {key}")

            return True

        except Exception as e:
            raise ValueError(f"Reference data integrity validation failed: {str(e)}")

    def get_dataset_info(self) -> Dict[str, Any]:
        """Get comprehensive information about the dataset."""
        metadata = self.load_metadata()
        index = self.load_index()

        # Count cases by category
        case_counts = {}
        for case_id, case_info in index["test_case_index"].items():
            for key, value in case_info.items():
                if key not in case_counts:
                    case_counts[key] = {}
                if value not in case_counts[key]:
                    case_counts[key][value] = 0
                case_counts[key][value] += 1

        return {
            "metadata": metadata,
            "total_cases": index["total_cases"],
            "case_distribution": case_counts,
            "reference_dir": str(self.reference_dir),
            "files_exist": {
                "reference_data": (self.reference_dir / "reference_data_v1.0.json.zst").exists(),
                "metadata": (self.reference_dir / "reference_metadata.json").exists(),
                "index": (self.reference_dir / "reference_index.json").exists()
            }
        }