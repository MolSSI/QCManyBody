"""
Validation reporting classes for ultra-strict numerical comparison.
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple, Union
import json


@dataclass
class ValidationError:
    """Details of a specific validation failure."""
    location: str
    expected_value: Any
    actual_value: Any
    absolute_difference: Optional[float] = None
    relative_difference: Optional[float] = None
    error_type: str = "numerical_mismatch"
    description: str = ""


@dataclass
class ValidationResult:
    """Result of a single validation check."""
    passed: bool
    error_message: str = ""
    absolute_difference: Optional[float] = None
    relative_difference: Optional[float] = None
    failure_indices: Optional[Tuple] = None
    max_absolute_difference: Optional[float] = None
    max_relative_difference: Optional[float] = None


@dataclass
class ValidationReport:
    """Comprehensive validation report for a test case."""
    test_id: str
    reference_key: str
    passed: bool
    timestamp: datetime = field(default_factory=datetime.utcnow)
    numerical_differences: Dict[str, float] = field(default_factory=dict)
    error_details: List[ValidationError] = field(default_factory=list)
    performance_metrics: Dict[str, float] = field(default_factory=dict)
    tolerance_used: float = 1e-12

    def add_error(self, location: str, expected: Any, actual: Any,
                  abs_diff: Optional[float] = None, rel_diff: Optional[float] = None,
                  error_type: str = "numerical_mismatch", description: str = "") -> None:
        """Add a validation error to the report."""
        error = ValidationError(
            location=location,
            expected_value=expected,
            actual_value=actual,
            absolute_difference=abs_diff,
            relative_difference=rel_diff,
            error_type=error_type,
            description=description
        )
        self.error_details.append(error)
        self.passed = False

    def add_numerical_difference(self, key: str, difference: float) -> None:
        """Add a numerical difference measurement."""
        self.numerical_differences[key] = difference

    def add_performance_metric(self, key: str, value: float) -> None:
        """Add a performance metric measurement."""
        self.performance_metrics[key] = value

    def generate_summary(self) -> str:
        """Generate a concise summary of the validation results."""
        if self.passed:
            return f"✅ PASS: {self.test_id} - All validations passed within {self.tolerance_used} tolerance"

        error_count = len(self.error_details)
        max_diff = max(self.numerical_differences.values()) if self.numerical_differences else 0

        return (f"❌ FAIL: {self.test_id} - {error_count} errors, "
                f"max difference: {max_diff:.2e}, tolerance: {self.tolerance_used}")

    def generate_detailed_report(self) -> str:
        """Generate a detailed validation report."""
        lines = [
            "=" * 80,
            f"VALIDATION REPORT: {self.test_id}",
            "=" * 80,
            f"Reference Key: {self.reference_key}",
            f"Timestamp: {self.timestamp.isoformat()}",
            f"Overall Result: {'PASS' if self.passed else 'FAIL'}",
            f"Tolerance Used: {self.tolerance_used}",
            "",
        ]

        if self.numerical_differences:
            lines.extend([
                "NUMERICAL DIFFERENCES:",
                "-" * 40,
            ])
            for key, diff in self.numerical_differences.items():
                status = "✅ OK" if diff <= self.tolerance_used else "❌ FAIL"
                lines.append(f"{status} {key}: {diff:.2e}")
            lines.append("")

        if self.error_details:
            lines.extend([
                f"VALIDATION ERRORS ({len(self.error_details)}):",
                "-" * 40,
            ])
            for i, error in enumerate(self.error_details, 1):
                lines.extend([
                    f"Error {i}: {error.error_type}",
                    f"  Location: {error.location}",
                    f"  Expected: {error.expected_value}",
                    f"  Actual: {error.actual_value}",
                ])
                if error.absolute_difference is not None:
                    lines.append(f"  Abs Diff: {error.absolute_difference:.2e}")
                if error.relative_difference is not None:
                    lines.append(f"  Rel Diff: {error.relative_difference:.2e}")
                if error.description:
                    lines.append(f"  Description: {error.description}")
                lines.append("")

        if self.performance_metrics:
            lines.extend([
                "PERFORMANCE METRICS:",
                "-" * 40,
            ])
            for key, value in self.performance_metrics.items():
                lines.append(f"{key}: {value:.4f}")
            lines.append("")

        lines.append("=" * 80)
        return "\n".join(lines)

    def save_to_file(self, filepath: str) -> None:
        """Save the validation report to a file."""
        report_data = {
            "test_id": self.test_id,
            "reference_key": self.reference_key,
            "passed": self.passed,
            "timestamp": self.timestamp.isoformat(),
            "tolerance_used": self.tolerance_used,
            "numerical_differences": self.numerical_differences,
            "performance_metrics": self.performance_metrics,
            "error_details": [
                {
                    "location": err.location,
                    "expected_value": str(err.expected_value),
                    "actual_value": str(err.actual_value),
                    "absolute_difference": err.absolute_difference,
                    "relative_difference": err.relative_difference,
                    "error_type": err.error_type,
                    "description": err.description
                }
                for err in self.error_details
            ]
        }

        with open(filepath, 'w') as f:
            json.dump(report_data, f, indent=2)


@dataclass
class BatchValidationReport:
    """Report for batch validation of multiple test cases."""
    batch_id: str
    total_cases: int
    passed_cases: int
    failed_cases: int
    timestamp: datetime = field(default_factory=datetime.utcnow)
    individual_reports: List[ValidationReport] = field(default_factory=list)
    overall_metrics: Dict[str, float] = field(default_factory=dict)

    @property
    def success_rate(self) -> float:
        """Calculate the success rate for the batch."""
        return self.passed_cases / self.total_cases if self.total_cases > 0 else 0.0

    def add_report(self, report: ValidationReport) -> None:
        """Add an individual validation report to the batch."""
        self.individual_reports.append(report)
        if report.passed:
            self.passed_cases += 1
        else:
            self.failed_cases += 1

    def generate_summary(self) -> str:
        """Generate a summary of the batch validation results."""
        return (f"Batch {self.batch_id}: {self.passed_cases}/{self.total_cases} passed "
                f"({self.success_rate:.1%} success rate)")

    def generate_detailed_report(self) -> str:
        """Generate a detailed batch validation report."""
        lines = [
            "=" * 100,
            f"BATCH VALIDATION REPORT: {self.batch_id}",
            "=" * 100,
            f"Timestamp: {self.timestamp.isoformat()}",
            f"Total Cases: {self.total_cases}",
            f"Passed: {self.passed_cases}",
            f"Failed: {self.failed_cases}",
            f"Success Rate: {self.success_rate:.1%}",
            "",
        ]

        if self.overall_metrics:
            lines.extend([
                "OVERALL METRICS:",
                "-" * 50,
            ])
            for key, value in self.overall_metrics.items():
                lines.append(f"{key}: {value:.4f}")
            lines.append("")

        # Summary of failed cases
        failed_reports = [r for r in self.individual_reports if not r.passed]
        if failed_reports:
            lines.extend([
                f"FAILED CASES ({len(failed_reports)}):",
                "-" * 50,
            ])
            for report in failed_reports:
                lines.append(f"❌ {report.test_id}: {len(report.error_details)} errors")

        lines.append("=" * 100)
        return "\n".join(lines)