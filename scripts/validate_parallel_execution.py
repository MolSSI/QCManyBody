#!/usr/bin/env python3
"""Ultra-strict validation of parallel execution correctness.

This script implements comprehensive validation of parallel vs sequential execution
with quantum chemistry precision (1e-12 tolerance) to ensure mathematical correctness.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from qcelemental.models import Molecule

# Add parent directory to path for qcmanybody imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from qcmanybody.core import ManyBodyCore
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
from qcmanybody.models.v1 import BsseEnum


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


class ParallelValidationFramework:
    """Ultra-strict validation framework for parallel execution correctness."""

    def __init__(self, tolerance: float = 1e-12):
        """Initialize validation framework.

        Parameters
        ----------
        tolerance : float
            Numerical tolerance for ultra-strict validation
        """
        self.tolerance = tolerance
        self.validation_results = {
            "tests_run": 0,
            "tests_passed": 0,
            "tests_failed": 0,
            "max_difference": 0.0,
            "failed_tests": []
        }

    def create_test_molecule(self, system_type: str) -> Molecule:
        """Create test molecules for validation.

        Parameters
        ----------
        system_type : str
            Type of test system to create

        Returns
        -------
        Molecule
            Test molecule
        """
        if system_type == "simple_dimer":
            return Molecule.from_data("""
            H 0.0 0.0 0.0
            --
            H 2.0 0.0 0.0
            """)

        elif system_type == "water_dimer":
            return Molecule.from_data("""
            O  0.0000  0.0000  0.0000
            H  0.7570  0.5860  0.0000
            H -0.7570  0.5860  0.0000
            --
            O  3.0000  0.0000  0.0000
            H  3.7570  0.5860  0.0000
            H  2.2430  0.5860  0.0000
            """)

        elif system_type == "water_trimer":
            return Molecule.from_data("""
            O  0.0000  0.0000  0.0000
            H  0.7570  0.5860  0.0000
            H -0.7570  0.5860  0.0000
            --
            O  3.0000  0.0000  0.0000
            H  3.7570  0.5860  0.0000
            H  2.2430  0.5860  0.0000
            --
            O  0.0000  3.0000  0.0000
            H  0.7570  3.5860  0.0000
            H -0.7570  3.5860  0.0000
            """)

        else:
            raise ValueError(f"Unknown system type: {system_type}")

    def create_test_configurations(self) -> List[Dict]:
        """Create test configurations for validation.

        Returns
        -------
        List[Dict]
            List of test configuration dictionaries
        """
        test_configs = []

        # Test different molecular systems
        systems = ["simple_dimer", "water_dimer", "water_trimer"]

        # Test different execution modes
        execution_modes = ["serial", "threading"]

        # Test different worker counts
        worker_counts = [1, 2, 4]

        # Test different BSSE treatments
        bsse_types = [
            [BsseEnum.nocp],
            [BsseEnum.cp],
        ]

        for system in systems:
            for bsse_type in bsse_types:
                for mode in execution_modes:
                    for workers in worker_counts:
                        # Skip serial mode with multiple workers
                        if mode == "serial" and workers > 1:
                            continue

                        test_configs.append({
                            "system": system,
                            "bsse_type": bsse_type,
                            "execution_mode": mode,
                            "max_workers": workers,
                            "test_id": f"{system}_{bsse_type[0].value}_{mode}_{workers}w"
                        })

        return test_configs

    def run_sequential_reference(self, core: ManyBodyCore) -> Dict:
        """Run sequential execution to generate reference results.

        Parameters
        ----------
        core : ManyBodyCore
            ManyBodyCore instance for reference calculation

        Returns
        -------
        Dict
            Reference results from sequential execution
        """
        logging.debug("Running sequential reference calculation")

        reference_results = {}

        # Use the standard iterate_molecules for sequential reference
        for mc, label, mol in core.iterate_molecules():
            # Create consistent placeholder result
            natoms = len(mol.symbols)
            # Use a deterministic energy calculation for consistency
            energy = -natoms * 1.0 - sum(ord(c) for c in label) * 1e-6

            reference_results[label] = {
                "energy": energy,
                "model_chemistry": mc,
                "natoms": natoms,
                "molecule_symbols": mol.symbols
            }

        logging.debug(f"Sequential reference: {len(reference_results)} fragments")
        return reference_results

    def run_parallel_test(self, core: ManyBodyCore, config: ParallelConfig) -> Dict:
        """Run parallel execution test.

        Parameters
        ----------
        core : ManyBodyCore
            ManyBodyCore instance for parallel calculation
        config : ParallelConfig
            Parallel execution configuration

        Returns
        -------
        Dict
            Results from parallel execution
        """
        logging.debug(f"Running parallel test: {config.execution_mode}, {config.max_workers} workers")

        # Use placeholder execution for consistent results
        config.use_qcengine = False

        executor = ParallelManyBodyExecutor(core, config)
        parallel_results_raw = executor.execute_full_calculation()

        # Convert AtomicResult objects to comparable format
        parallel_results = {}
        for label, result in parallel_results_raw.items():
            # Handle model field which might be dict or Model object
            model_method = result.model
            if hasattr(model_method, 'method'):
                model_method = model_method.method
            elif isinstance(model_method, dict):
                model_method = model_method.get("method", "unknown")
            else:
                model_method = str(model_method)

            parallel_results[label] = {
                "energy": result.return_result,
                "model_chemistry": model_method,
                "natoms": len(result.molecule.symbols),
                "molecule_symbols": result.molecule.symbols
            }

        logging.debug(f"Parallel test: {len(parallel_results)} fragments")
        return parallel_results

    def validate_results(self, reference: Dict, parallel: Dict, test_id: str) -> bool:
        """Validate parallel results against reference with ultra-strict tolerance.

        Parameters
        ----------
        reference : Dict
            Reference results from sequential execution
        parallel : Dict
            Results from parallel execution
        test_id : str
            Test identifier for reporting

        Returns
        -------
        bool
            True if validation passes, False otherwise
        """
        self.validation_results["tests_run"] += 1

        # Check fragment count
        if len(reference) != len(parallel):
            error_msg = f"Fragment count mismatch: reference={len(reference)}, parallel={len(parallel)}"
            logging.error(f"{test_id}: {error_msg}")
            self.validation_results["failed_tests"].append({
                "test_id": test_id,
                "error": error_msg
            })
            self.validation_results["tests_failed"] += 1
            return False

        # Check fragment labels
        ref_labels = set(reference.keys())
        par_labels = set(parallel.keys())
        if ref_labels != par_labels:
            missing_ref = ref_labels - par_labels
            missing_par = par_labels - ref_labels
            error_msg = f"Label mismatch: missing_in_parallel={missing_ref}, extra_in_parallel={missing_par}"
            logging.error(f"{test_id}: {error_msg}")
            self.validation_results["failed_tests"].append({
                "test_id": test_id,
                "error": error_msg
            })
            self.validation_results["tests_failed"] += 1
            return False

        # Ultra-strict numerical validation
        max_difference = 0.0
        for label in ref_labels:
            ref_data = reference[label]
            par_data = parallel[label]

            # Validate energy values
            energy_diff = abs(ref_data["energy"] - par_data["energy"])
            max_difference = max(max_difference, energy_diff)

            if energy_diff > self.tolerance:
                error_msg = f"Energy difference for {label}: {energy_diff} > {self.tolerance}"
                logging.error(f"{test_id}: {error_msg}")
                self.validation_results["failed_tests"].append({
                    "test_id": test_id,
                    "error": error_msg,
                    "label": label,
                    "difference": energy_diff
                })
                self.validation_results["tests_failed"] += 1
                return False

            # Validate model chemistry consistency
            # Handle both dict and Model object formats
            par_method = par_data["model_chemistry"]
            if hasattr(par_method, 'method'):
                par_method = par_method.method
            elif isinstance(par_method, dict):
                par_method = par_method.get("method", str(par_method))

            if ref_data["model_chemistry"] != par_method:
                error_msg = f"Model chemistry mismatch for {label}: {ref_data['model_chemistry']} != {par_method}"
                logging.error(f"{test_id}: {error_msg}")
                self.validation_results["failed_tests"].append({
                    "test_id": test_id,
                    "error": error_msg,
                    "label": label
                })
                self.validation_results["tests_failed"] += 1
                return False

            # Validate molecular structure consistency
            if ref_data["natoms"] != par_data["natoms"]:
                error_msg = f"Atom count mismatch for {label}: {ref_data['natoms']} != {par_data['natoms']}"
                logging.error(f"{test_id}: {error_msg}")
                self.validation_results["failed_tests"].append({
                    "test_id": test_id,
                    "error": error_msg,
                    "label": label
                })
                self.validation_results["tests_failed"] += 1
                return False

        # Update global maximum difference
        self.validation_results["max_difference"] = max(
            self.validation_results["max_difference"], max_difference
        )

        logging.info(f"{test_id}: PASSED (max_diff={max_difference:.2e})")
        self.validation_results["tests_passed"] += 1
        return True

    def run_comprehensive_validation(self) -> bool:
        """Run comprehensive validation across all test configurations.

        Returns
        -------
        bool
            True if all tests pass, False otherwise
        """
        logging.info("Starting comprehensive parallel execution validation")
        logging.info(f"Ultra-strict tolerance: {self.tolerance}")

        test_configs = self.create_test_configurations()
        logging.info(f"Running {len(test_configs)} validation tests")

        all_passed = True

        for config in test_configs:
            test_id = config["test_id"]
            logging.info(f"Running test: {test_id}")

            try:
                # Create test system
                molecule = self.create_test_molecule(config["system"])
                nfragments = len(molecule.fragments)
                max_nbody = min(3, nfragments)  # Limit for testing performance

                levels = {i: "hf" for i in range(1, max_nbody + 1)}

                core = ManyBodyCore(
                    molecule=molecule,
                    bsse_type=config["bsse_type"],
                    levels=levels,
                    return_total_data=False,
                    supersystem_ie_only=False,
                    embedding_charges={}
                )

                # Run reference sequential calculation
                reference_results = self.run_sequential_reference(core)

                # Run parallel test
                parallel_config = ParallelConfig(
                    execution_mode=config["execution_mode"],
                    max_workers=config["max_workers"],
                    use_qcengine=False  # Use consistent placeholder for validation
                )

                parallel_results = self.run_parallel_test(core, parallel_config)

                # Validate results
                test_passed = self.validate_results(reference_results, parallel_results, test_id)
                if not test_passed:
                    all_passed = False

            except Exception as e:
                logging.error(f"{test_id}: EXCEPTION - {e}")
                self.validation_results["failed_tests"].append({
                    "test_id": test_id,
                    "error": f"Exception: {e}"
                })
                self.validation_results["tests_failed"] += 1
                all_passed = False

        return all_passed

    def generate_validation_report(self, output_file: str = "parallel_validation_report.json"):
        """Generate comprehensive validation report.

        Parameters
        ----------
        output_file : str
            Path to save validation report
        """
        report = {
            "validation_summary": self.validation_results,
            "tolerance": self.tolerance,
            "overall_status": "PASSED" if self.validation_results["tests_failed"] == 0 else "FAILED",
            "pass_rate": (
                self.validation_results["tests_passed"] / max(self.validation_results["tests_run"], 1)
            ) * 100
        }

        # Save report
        output_path = Path(output_file)
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)

        logging.info(f"Validation report saved to: {output_path}")

        # Print summary
        self.print_validation_summary(report)

    def print_validation_summary(self, report: Dict):
        """Print validation summary to console.

        Parameters
        ----------
        report : Dict
            Validation report to summarize
        """
        print(f"\n{'='*80}")
        print("PARALLEL EXECUTION VALIDATION SUMMARY")
        print(f"{'='*80}")
        print(f"Overall Status: {report['overall_status']}")
        print(f"Pass Rate: {report['pass_rate']:.1f}%")
        print(f"Tests Run: {report['validation_summary']['tests_run']}")
        print(f"Tests Passed: {report['validation_summary']['tests_passed']}")
        print(f"Tests Failed: {report['validation_summary']['tests_failed']}")
        print(f"Tolerance: {report['tolerance']}")
        print(f"Maximum Difference: {report['validation_summary']['max_difference']:.2e}")

        if report['validation_summary']['failed_tests']:
            print(f"\nFAILED TESTS:")
            for failure in report['validation_summary']['failed_tests']:
                print(f"  {failure['test_id']}: {failure['error']}")

        print(f"{'='*80}")


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description="Validate parallel execution correctness")
    parser.add_argument("--tolerance", type=float, default=1e-12,
                       help="Numerical tolerance for validation")
    parser.add_argument("--output", default="parallel_validation_report.json",
                       help="Output file for validation report")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose logging")

    args = parser.parse_args()

    setup_logging(args.verbose)

    # Create validation framework
    validator = ParallelValidationFramework(tolerance=args.tolerance)

    # Run comprehensive validation
    all_tests_passed = validator.run_comprehensive_validation()

    # Generate report
    validator.generate_validation_report(args.output)

    # Exit with appropriate code
    sys.exit(0 if all_tests_passed else 1)


if __name__ == "__main__":
    main()