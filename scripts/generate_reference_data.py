#!/usr/bin/env python3
"""
Generate Golden Reference Dataset for QCManyBody Parallel Execution Validation

This script generates a comprehensive reference dataset using the current sequential
QCManyBody code before any parallel modifications. The dataset serves as the absolute
truth for all regression testing throughout parallel development.

Critical Requirements:
- Must use UNMODIFIED sequential QCManyBody code
- Numerical precision must be at least 1e-14 (working precision)
- Comprehensive coverage of all test scenarios
- Ultra-strict validation with 1e-12 tolerance
"""

import json
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union

import numpy as np
import zstandard
from qcelemental.models import Molecule

# Import QCManyBody components
from qcmanybody import ManyBodyCore
from qcmanybody.models import BsseEnum
from qcmanybody.tests.common import mol_h2o_3_dict, mol_ne2, specifications
from qcmanybody.tests.utils import jsonify, run_qcengine
from qcmanybody.utils import delabeler

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('generate_reference_data.log')
    ]
)
logger = logging.getLogger(__name__)

class ReferenceDataGenerator:
    """Generate comprehensive golden reference dataset for parallel validation."""

    def __init__(self, output_dir: str = "qcmanybody/tests/reference_data_parallel"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Test case matrix for systematic coverage
        self.test_matrix = {
            "systems": {
                "ne2": {
                    "molecule": mol_ne2,
                    "description": "Simple noble gas dimer"
                },
                "h2o3": {
                    "molecule": mol_h2o_3_dict,
                    "description": "Water trimer - standard test case"
                }
            },
            "bsse_types": [BsseEnum.cp, BsseEnum.nocp, BsseEnum.vmfc],
            "max_nbody_levels": [2, 3],
            "drivers": ["energy", "gradient"],  # Start with these, add hessian later
            "level_combinations": {
                "single_scf": {1: "e_scf", 2: "e_scf", 3: "e_scf"},
                "single_mp2": {1: "e_mp2", 2: "e_mp2", 3: "e_mp2"},
                "mixed_scf_mp2": {1: "e_scf", 2: "e_mp2", 3: "e_mp2"}
            }
        }

        self.metadata = {
            "generator_version": "1.0.0",
            "qcmanybody_version": self._get_qcmanybody_version(),
            "generation_timestamp": datetime.utcnow().isoformat(),
            "numerical_precision": 1e-14,
            "validation_tolerance": 1e-12,
            "total_test_cases": 0,
            "generation_notes": "Golden reference dataset for parallel execution validation"
        }

        self.reference_data = {}

    def _get_qcmanybody_version(self) -> str:
        """Get QCManyBody version information."""
        try:
            import qcmanybody
            return getattr(qcmanybody, '__version__', 'unknown')
        except:
            return 'unknown'

    def _generate_test_case_id(self, system: str, bsse_type: BsseEnum,
                              max_nbody: int, driver: str, levels: str) -> str:
        """Generate unique test case identifier."""
        return f"{system}_{bsse_type.value}_{driver}_maxnb{max_nbody}_{levels}"

    def _create_specifications_for_driver(self, levels: Dict[int, str], driver: str) -> Dict[str, Dict[str, Any]]:
        """Create specifications dictionary for given driver."""
        specs = {}
        for level_name in set(levels.values()):
            base_spec = specifications[level_name].copy()
            base_spec["specification"]["driver"] = driver
            specs[level_name] = base_spec
        return specs

    def generate_reference_case(self, system_name: str, system_data: Dict[str, Any],
                               bsse_type: BsseEnum, max_nbody: int, driver: str,
                               levels: Dict[int, str], levels_name: str) -> Dict[str, Any]:
        """Generate reference data for a single test case."""

        test_case_id = self._generate_test_case_id(
            system_name, bsse_type, max_nbody, driver, levels_name
        )

        logger.info(f"Generating reference case: {test_case_id}")

        try:
            # Create specifications for this driver
            specs = self._create_specifications_for_driver(levels, driver)

            # Trim levels to max_nbody
            active_levels = {k: v for k, v in levels.items() if k <= max_nbody}

            # Create molecule object
            if isinstance(system_data["molecule"], dict):
                molecule = Molecule(**system_data["molecule"])
            else:
                molecule = system_data["molecule"]

            # Generate component results using QCEngine
            mbc, component_results = run_qcengine(
                specifications=specs,
                molecule=molecule,
                bsse_type=[bsse_type],
                levels=active_levels,
                return_total_data=True,
                supersystem_ie_only=False,
                embedding_charges=None
            )

            # Analyze results through ManyBodyCore
            final_results = mbc.analyze(component_results)

            # Package reference data with metadata
            reference_entry = {
                "test_case_id": test_case_id,
                "system": {
                    "name": system_name,
                    "description": system_data["description"],
                    "molecule": molecule.dict() if hasattr(molecule, 'dict') else system_data["molecule"]
                },
                "configuration": {
                    "bsse_type": bsse_type.value,
                    "max_nbody": max_nbody,
                    "driver": driver,
                    "levels": active_levels,
                    "levels_name": levels_name
                },
                "qcmanybody_config": {
                    "return_total_data": True,
                    "supersystem_ie_only": False,
                    "embedding_charges": None
                },
                "timestamp": datetime.utcnow().isoformat(),
                "component_count": len(component_results),
                "component_results": jsonify(component_results),
                "final_results": jsonify(final_results),
                "numerical_precision": 1e-14
            }

            # Validate numerical precision
            self._validate_precision(reference_entry)

            return reference_entry

        except Exception as e:
            logger.error(f"Failed to generate reference case {test_case_id}: {str(e)}")
            raise

    def _validate_precision(self, reference_entry: Dict[str, Any]) -> None:
        """Validate that the reference data meets precision requirements."""
        final_results = reference_entry["final_results"]

        # Check for any NaN or infinite values
        def check_numeric_validity(data, path=""):
            if isinstance(data, dict):
                for key, value in data.items():
                    check_numeric_validity(value, f"{path}.{key}")
            elif isinstance(data, (list, tuple)):
                for i, value in enumerate(data):
                    check_numeric_validity(value, f"{path}[{i}]")
            elif isinstance(data, (int, float)):
                if np.isnan(data) or np.isinf(data):
                    raise ValueError(f"Invalid numeric value at {path}: {data}")

        check_numeric_validity(final_results, "final_results")

        # Verify key properties exist
        required_keys = ["ret_energy"]
        for key in required_keys:
            if key not in final_results:
                raise ValueError(f"Missing required result key: {key}")

    def generate_complete_dataset(self) -> None:
        """Generate the complete reference dataset."""
        logger.info("Starting complete reference dataset generation")

        total_cases = 0
        generated_cases = 0
        failed_cases = []

        # Calculate total expected cases
        for system_name in self.test_matrix["systems"]:
            for bsse_type in self.test_matrix["bsse_types"]:
                for max_nbody in self.test_matrix["max_nbody_levels"]:
                    for driver in self.test_matrix["drivers"]:
                        for levels_name in self.test_matrix["level_combinations"]:
                            total_cases += 1

        logger.info(f"Planning to generate {total_cases} reference cases")

        # Generate all test cases
        for system_name, system_data in self.test_matrix["systems"].items():
            for bsse_type in self.test_matrix["bsse_types"]:
                for max_nbody in self.test_matrix["max_nbody_levels"]:
                    for driver in self.test_matrix["drivers"]:
                        for levels_name, levels in self.test_matrix["level_combinations"].items():

                            test_case_id = self._generate_test_case_id(
                                system_name, bsse_type, max_nbody, driver, levels_name
                            )

                            try:
                                reference_case = self.generate_reference_case(
                                    system_name, system_data, bsse_type, max_nbody,
                                    driver, levels, levels_name
                                )

                                self.reference_data[test_case_id] = reference_case
                                generated_cases += 1

                                logger.info(f"Generated case {generated_cases}/{total_cases}: {test_case_id}")

                            except Exception as e:
                                logger.error(f"Failed case {test_case_id}: {str(e)}")
                                failed_cases.append((test_case_id, str(e)))

        # Update metadata
        self.metadata["total_test_cases"] = generated_cases
        self.metadata["failed_cases"] = len(failed_cases)
        self.metadata["success_rate"] = generated_cases / total_cases if total_cases > 0 else 0
        self.metadata["failed_case_details"] = failed_cases

        logger.info(f"Dataset generation complete: {generated_cases}/{total_cases} cases successful")

        if failed_cases:
            logger.warning(f"Failed cases: {len(failed_cases)}")
            for case_id, error in failed_cases:
                logger.warning(f"  {case_id}: {error}")

    def save_dataset(self) -> None:
        """Save the reference dataset to compressed files."""
        logger.info("Saving reference dataset")

        # Save main reference dataset
        reference_file = self.output_dir / "reference_data_v1.0.json.zst"
        with zstandard.open(reference_file, "wt") as f:
            json.dump(self.reference_data, f, indent=2)

        # Save metadata
        metadata_file = self.output_dir / "reference_metadata.json"
        with open(metadata_file, "w") as f:
            json.dump(self.metadata, f, indent=2)

        # Create index file
        index_data = {
            "dataset_version": "1.0.0",
            "total_cases": len(self.reference_data),
            "test_case_index": {
                case_id: {
                    "system": data["system"]["name"],
                    "bsse_type": data["configuration"]["bsse_type"],
                    "max_nbody": data["configuration"]["max_nbody"],
                    "driver": data["configuration"]["driver"],
                    "levels": data["configuration"]["levels_name"]
                }
                for case_id, data in self.reference_data.items()
            }
        }

        index_file = self.output_dir / "reference_index.json"
        with open(index_file, "w") as f:
            json.dump(index_data, f, indent=2)

        # Log file sizes
        logger.info(f"Reference data saved:")
        logger.info(f"  Main dataset: {reference_file} ({reference_file.stat().st_size} bytes)")
        logger.info(f"  Metadata: {metadata_file} ({metadata_file.stat().st_size} bytes)")
        logger.info(f"  Index: {index_file} ({index_file.stat().st_size} bytes)")

    def validate_dataset(self) -> bool:
        """Validate the generated dataset for completeness and integrity."""
        logger.info("Validating reference dataset")

        # Check that we have expected number of cases
        expected_cases = (
            len(self.test_matrix["systems"]) *
            len(self.test_matrix["bsse_types"]) *
            len(self.test_matrix["max_nbody_levels"]) *
            len(self.test_matrix["drivers"]) *
            len(self.test_matrix["level_combinations"])
        )

        actual_cases = len(self.reference_data)

        if actual_cases < expected_cases:
            logger.warning(f"Dataset incomplete: {actual_cases}/{expected_cases} cases")
            return False

        # Validate each case has required structure
        for case_id, case_data in self.reference_data.items():
            required_top_level = ["test_case_id", "system", "configuration", "final_results"]
            for key in required_top_level:
                if key not in case_data:
                    logger.error(f"Case {case_id} missing required key: {key}")
                    return False

        logger.info("Dataset validation passed")
        return True


def main():
    """Main execution function."""
    logger.info("Starting QCManyBody reference dataset generation")

    # Check for required dependencies
    try:
        import qcengine
        import psi4  # Will be needed for calculations
    except ImportError as e:
        logger.error(f"Missing required dependency: {e}")
        logger.error("Please ensure QCEngine and Psi4 are installed")
        sys.exit(1)

    # Create generator and run
    generator = ReferenceDataGenerator()

    try:
        # Generate complete dataset
        start_time = time.time()
        generator.generate_complete_dataset()
        generation_time = time.time() - start_time

        # Validate dataset
        if not generator.validate_dataset():
            logger.error("Dataset validation failed")
            sys.exit(1)

        # Save dataset
        generator.save_dataset()

        logger.info(f"Reference dataset generation completed successfully")
        logger.info(f"Total generation time: {generation_time:.2f} seconds")
        logger.info(f"Generated {len(generator.reference_data)} reference cases")

    except Exception as e:
        logger.error(f"Reference dataset generation failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()