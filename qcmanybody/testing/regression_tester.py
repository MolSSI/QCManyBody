"""
Ultra-strict regression testing framework for QCManyBody parallel execution validation.

Implements 1e-12 tolerance numerical comparison for ensuring parallel execution
produces identical results to sequential execution.
"""

import math
import time
from typing import Any, Dict, List, Optional, Union

import numpy as np

from .reference_loader import ReferenceDataLoader
from .validation_report import ValidationReport, ValidationResult, BatchValidationReport


class ParallelRegressionTester:
    """Ultra-strict regression tester for parallel execution validation."""

    def __init__(self, tolerance: float = 1e-12, reference_data_path: Optional[str] = None):
        """Initialize the regression tester.

        Args:
            tolerance: Numerical tolerance for comparisons (default 1e-12)
            reference_data_path: Path to reference data directory
        """
        self.tolerance = tolerance
        self.reference_loader = ReferenceDataLoader(reference_data_path)

    def validate_result(self, test_result: Dict[str, Any], reference_key: str) -> ValidationReport:
        """Validate a test result against reference data.

        Args:
            test_result: Result from parallel execution to validate
            reference_key: Key to identify reference case

        Returns:
            ValidationReport with detailed comparison results
        """
        start_time = time.time()

        # Create validation report
        report = ValidationReport(
            test_id=reference_key,
            reference_key=reference_key,
            passed=True,
            tolerance_used=self.tolerance
        )

        try:
            # Load reference data
            reference_case = self.reference_loader.get_reference_case(reference_key)
            reference_results = reference_case["final_results"]

            # Validate core results
            self._validate_final_results(test_result, reference_results, report)

        except Exception as e:
            report.passed = False
            report.add_error(
                location="validation_process",
                expected="successful_validation",
                actual=f"exception_{type(e).__name__}",
                error_type="validation_error",
                description=str(e)
            )

        # Add performance metrics
        validation_time = time.time() - start_time
        report.add_performance_metric("validation_time_seconds", validation_time)

        return report

    def validate_batch(self, test_results: Dict[str, Dict[str, Any]],
                      reference_keys: List[str]) -> BatchValidationReport:
        """Validate multiple test results in batch.

        Args:
            test_results: Dictionary of test results keyed by test case ID
            reference_keys: List of reference keys to validate against

        Returns:
            BatchValidationReport with summary and individual results
        """
        batch_report = BatchValidationReport(
            batch_id=f"batch_{int(time.time())}",
            total_cases=len(reference_keys),
            passed_cases=0,
            failed_cases=0
        )

        start_time = time.time()

        for ref_key in reference_keys:
            if ref_key in test_results:
                individual_report = self.validate_result(test_results[ref_key], ref_key)
                batch_report.add_report(individual_report)
            else:
                # Missing test result
                missing_report = ValidationReport(
                    test_id=ref_key,
                    reference_key=ref_key,
                    passed=False,
                    tolerance_used=self.tolerance
                )
                missing_report.add_error(
                    location="test_execution",
                    expected="test_result_present",
                    actual="test_result_missing",
                    error_type="missing_result",
                    description=f"No test result provided for {ref_key}"
                )
                batch_report.add_report(missing_report)

        # Add overall metrics
        batch_time = time.time() - start_time
        batch_report.overall_metrics["total_validation_time"] = batch_time
        batch_report.overall_metrics["average_time_per_case"] = batch_time / len(reference_keys)

        return batch_report

    def _validate_final_results(self, test_result: Dict[str, Any],
                               reference_result: Dict[str, Any],
                               report: ValidationReport) -> None:
        """Validate final results with ultra-strict precision."""

        # Core energy validation
        if "ret_energy" in reference_result:
            if "ret_energy" not in test_result:
                report.add_error(
                    location="ret_energy",
                    expected="present",
                    actual="missing",
                    error_type="missing_property"
                )
            else:
                validation_result = self.compare_floats(
                    test_result["ret_energy"],
                    reference_result["ret_energy"]
                )
                if not validation_result.passed:
                    report.add_error(
                        location="ret_energy",
                        expected=reference_result["ret_energy"],
                        actual=test_result["ret_energy"],
                        abs_diff=validation_result.absolute_difference,
                        rel_diff=validation_result.relative_difference,
                        error_type="energy_mismatch"
                    )
                else:
                    report.add_numerical_difference(
                        "ret_energy_abs_diff",
                        validation_result.absolute_difference or 0.0
                    )

        # Gradient validation
        if "ret_gradient" in reference_result:
            if "ret_gradient" not in test_result:
                report.add_error(
                    location="ret_gradient",
                    expected="present",
                    actual="missing",
                    error_type="missing_property"
                )
            else:
                validation_result = self.compare_arrays(
                    np.array(test_result["ret_gradient"]),
                    np.array(reference_result["ret_gradient"])
                )
                if not validation_result.passed:
                    report.add_error(
                        location="ret_gradient",
                        expected="array_match",
                        actual="array_mismatch",
                        abs_diff=validation_result.max_absolute_difference,
                        rel_diff=validation_result.max_relative_difference,
                        error_type="gradient_mismatch",
                        description=f"Failed elements at indices: {validation_result.failure_indices}"
                    )
                else:
                    report.add_numerical_difference(
                        "ret_gradient_max_abs_diff",
                        validation_result.max_absolute_difference or 0.0
                    )

        # Hessian validation
        if "ret_hessian" in reference_result:
            if "ret_hessian" not in test_result:
                report.add_error(
                    location="ret_hessian",
                    expected="present",
                    actual="missing",
                    error_type="missing_property"
                )
            else:
                validation_result = self.compare_arrays(
                    np.array(test_result["ret_hessian"]),
                    np.array(reference_result["ret_hessian"])
                )
                if not validation_result.passed:
                    report.add_error(
                        location="ret_hessian",
                        expected="array_match",
                        actual="array_mismatch",
                        abs_diff=validation_result.max_absolute_difference,
                        rel_diff=validation_result.max_relative_difference,
                        error_type="hessian_mismatch"
                    )
                else:
                    report.add_numerical_difference(
                        "ret_hessian_max_abs_diff",
                        validation_result.max_absolute_difference or 0.0
                    )

        # Energy body dictionary validation
        if "energy_body_dict" in reference_result:
            self._validate_body_dict(
                test_result.get("energy_body_dict", {}),
                reference_result["energy_body_dict"],
                "energy_body_dict",
                report
            )

        # Gradient body dictionary validation
        if "gradient_body_dict" in reference_result:
            self._validate_body_dict(
                test_result.get("gradient_body_dict", {}),
                reference_result["gradient_body_dict"],
                "gradient_body_dict",
                report,
                is_array=True
            )

        # Results properties validation
        if "results" in reference_result:
            self._validate_results_dict(
                test_result.get("results", {}),
                reference_result["results"],
                report
            )

    def _validate_body_dict(self, test_dict: Dict[str, Any], ref_dict: Dict[str, Any],
                           dict_name: str, report: ValidationReport, is_array: bool = False) -> None:
        """Validate body dictionaries (energy_body_dict, gradient_body_dict)."""

        for key, ref_value in ref_dict.items():
            # Handle string/int key conversion (JSON serialization converts int keys to strings)
            test_key = key
            if key not in test_dict:
                # Try converting between string and int keys
                if isinstance(key, str) and key.isdigit():
                    test_key = int(key)
                elif isinstance(key, int):
                    test_key = str(key)

                if test_key not in test_dict:
                    report.add_error(
                        location=f"{dict_name}.{key}",
                        expected="present",
                        actual="missing",
                        error_type="missing_key"
                    )
                    continue

            test_value = test_dict[test_key]

            if isinstance(ref_value, dict):
                # Nested dictionary - recurse
                if not isinstance(test_value, dict):
                    report.add_error(
                        location=f"{dict_name}.{key}",
                        expected="dict",
                        actual=type(test_value).__name__,
                        error_type="type_mismatch"
                    )
                else:
                    self._validate_body_dict(test_value, ref_value, f"{dict_name}.{key}", report, is_array)

            elif is_array and isinstance(ref_value, (list, np.ndarray)):
                # Array comparison
                validation_result = self.compare_arrays(np.array(test_value), np.array(ref_value))
                if not validation_result.passed:
                    report.add_error(
                        location=f"{dict_name}.{key}",
                        expected="array_match",
                        actual="array_mismatch",
                        abs_diff=validation_result.max_absolute_difference,
                        rel_diff=validation_result.max_relative_difference,
                        error_type="array_mismatch"
                    )

            else:
                # Scalar comparison
                validation_result = self.compare_floats(test_value, ref_value)
                if not validation_result.passed:
                    report.add_error(
                        location=f"{dict_name}.{key}",
                        expected=ref_value,
                        actual=test_value,
                        abs_diff=validation_result.absolute_difference,
                        rel_diff=validation_result.relative_difference,
                        error_type="scalar_mismatch"
                    )

    def _validate_results_dict(self, test_results: Dict[str, Any],
                              ref_results: Dict[str, Any],
                              report: ValidationReport) -> None:
        """Validate the results dictionary containing QC variables."""

        for key, ref_value in ref_results.items():
            if key not in test_results:
                report.add_error(
                    location=f"results.{key}",
                    expected="present",
                    actual="missing",
                    error_type="missing_qc_variable"
                )
                continue

            test_value = test_results[key]

            if isinstance(ref_value, (list, np.ndarray)):
                validation_result = self.compare_arrays(np.array(test_value), np.array(ref_value))
                if not validation_result.passed:
                    report.add_error(
                        location=f"results.{key}",
                        expected="array_match",
                        actual="array_mismatch",
                        abs_diff=validation_result.max_absolute_difference,
                        rel_diff=validation_result.max_relative_difference,
                        error_type="qc_variable_array_mismatch"
                    )
            else:
                validation_result = self.compare_floats(test_value, ref_value)
                if not validation_result.passed:
                    report.add_error(
                        location=f"results.{key}",
                        expected=ref_value,
                        actual=test_value,
                        abs_diff=validation_result.absolute_difference,
                        rel_diff=validation_result.relative_difference,
                        error_type="qc_variable_mismatch"
                    )

    def compare_floats(self, a: float, b: float, tolerance: Optional[float] = None) -> ValidationResult:
        """Ultra-strict floating point comparison.

        Args:
            a: First value
            b: Second value
            tolerance: Override default tolerance

        Returns:
            ValidationResult with comparison details
        """
        if tolerance is None:
            tolerance = self.tolerance

        # Handle special cases
        if math.isnan(a) or math.isnan(b):
            return ValidationResult(
                passed=False,
                error_message="NaN values detected"
            )

        if math.isinf(a) or math.isinf(b):
            # Both infinite with same sign is OK
            if math.isinf(a) and math.isinf(b) and (a > 0) == (b > 0):
                return ValidationResult(passed=True)
            else:
                return ValidationResult(
                    passed=False,
                    error_message="Infinite values detected or sign mismatch"
                )

        # Calculate differences
        abs_diff = abs(a - b)

        # Avoid division by zero for relative difference
        denominator = max(abs(a), abs(b), 1e-15)
        rel_diff = abs_diff / denominator

        # Pass if both absolute and relative differences are within tolerance
        passed = (abs_diff <= tolerance) or (rel_diff <= tolerance)

        return ValidationResult(
            passed=passed,
            absolute_difference=abs_diff,
            relative_difference=rel_diff,
            error_message="" if passed else f"Numerical mismatch: {a} vs {b}"
        )

    def compare_arrays(self, a: np.ndarray, b: np.ndarray,
                      tolerance: Optional[float] = None) -> ValidationResult:
        """Ultra-strict array comparison with element-wise analysis.

        Args:
            a: First array
            b: Second array
            tolerance: Override default tolerance

        Returns:
            ValidationResult with detailed array comparison
        """
        if tolerance is None:
            tolerance = self.tolerance

        # Shape check
        if a.shape != b.shape:
            return ValidationResult(
                passed=False,
                error_message=f"Shape mismatch: {a.shape} vs {b.shape}"
            )

        # Convert to float arrays for numerical comparison
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)

        # Check for NaN or infinite values
        if np.any(np.isnan(a)) or np.any(np.isnan(b)):
            return ValidationResult(
                passed=False,
                error_message="NaN values detected in arrays"
            )

        if np.any(np.isinf(a)) or np.any(np.isinf(b)):
            # Check if infinite values match
            inf_a = np.isinf(a)
            inf_b = np.isinf(b)
            if not np.array_equal(inf_a, inf_b) or not np.array_equal(a[inf_a], b[inf_b]):
                return ValidationResult(
                    passed=False,
                    error_message="Infinite values don't match between arrays"
                )

        # Element-wise comparison
        abs_diff = np.abs(a - b)

        # Relative difference with safe division
        denominator = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-15)
        rel_diff = abs_diff / denominator

        # Elements fail if both absolute and relative differences exceed tolerance
        abs_failures = abs_diff > tolerance
        rel_failures = rel_diff > tolerance
        failures = abs_failures & rel_failures

        if np.any(failures):
            failure_indices = np.where(failures)
            max_abs_diff = np.max(abs_diff[failures])
            max_rel_diff = np.max(rel_diff[failures])

            return ValidationResult(
                passed=False,
                error_message=f"Array mismatch: {np.sum(failures)} failed elements",
                failure_indices=failure_indices,
                max_absolute_difference=max_abs_diff,
                max_relative_difference=max_rel_diff
            )

        # All elements passed
        max_abs_diff = np.max(abs_diff) if abs_diff.size > 0 else 0.0
        max_rel_diff = np.max(rel_diff) if rel_diff.size > 0 else 0.0

        return ValidationResult(
            passed=True,
            max_absolute_difference=max_abs_diff,
            max_relative_difference=max_rel_diff
        )

    def compare_energies(self, result_energy: float, reference_energy: float) -> bool:
        """Convenience method for energy comparison."""
        return self.compare_floats(result_energy, reference_energy).passed

    def compare_gradients(self, result_grad: np.ndarray, reference_grad: np.ndarray) -> bool:
        """Convenience method for gradient comparison."""
        return self.compare_arrays(result_grad, reference_grad).passed

    def compare_hessians(self, result_hess: np.ndarray, reference_hess: np.ndarray) -> bool:
        """Convenience method for hessian comparison."""
        return self.compare_arrays(result_hess, reference_hess).passed

    def compare_properties(self, result_props: Dict[str, Any],
                          reference_props: Dict[str, Any]) -> bool:
        """Convenience method for properties comparison."""
        report = ValidationReport(
            test_id="properties_comparison",
            reference_key="properties_comparison",
            passed=True,
            tolerance_used=self.tolerance
        )

        self._validate_results_dict(result_props, reference_props, report)
        return report.passed