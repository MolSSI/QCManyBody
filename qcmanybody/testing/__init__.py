"""
QCManyBody Testing Infrastructure

This module provides ultra-strict validation and regression testing infrastructure
for parallel execution development with 1e-12 numerical tolerance.
"""

from .regression_tester import ParallelRegressionTester
from .reference_loader import ReferenceDataLoader
from .validation_report import ValidationReport, ValidationError

__all__ = [
    "ParallelRegressionTester",
    "ReferenceDataLoader",
    "ValidationReport",
    "ValidationError"
]