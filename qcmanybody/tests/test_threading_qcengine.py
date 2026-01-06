"""Test threading-based QCEngine integration.

This module tests that QCEngine works correctly in multithreaded environments,
addressing the race condition issues discovered in ThreadPoolExecutor contexts.
"""

import pytest
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple

import qcengine as qcng
import qcelemental as qcel
from qcelemental.models import AtomicInput, AtomicResult, Molecule

# Skip all tests if QCEngine or Psi4 not available
HAS_QCENGINE = True
try:
    import qcengine as qcng
except ImportError:
    HAS_QCENGINE = False

HAS_PSI4 = False
if HAS_QCENGINE:
    try:
        # Test if Psi4 is available through QCEngine
        test_mol = qcel.models.Molecule.from_data("H 0 0 0\nH 0 0 1")
        test_input = qcel.models.AtomicInput(
            molecule=test_mol,
            driver="energy",
            model={"method": "hf", "basis": "sto-3g"}
        )
        test_result = qcng.compute(test_input, "psi4", raise_error=False)
        HAS_PSI4 = test_result.success
    except Exception:
        HAS_PSI4 = False

pytestmark = pytest.mark.skipif(
    not HAS_QCENGINE or not HAS_PSI4,
    reason="QCEngine and Psi4 required for threading integration tests"
)


class ThreadingQCEngineValidator:
    """Validator for QCEngine threading compatibility."""

    @staticmethod
    def create_test_molecule() -> Molecule:
        """Create simple water molecule for testing."""
        return qcel.models.Molecule.from_data("""
        O 0.0000000000  0.0000000000 -0.0657755706
        H 0.0000000000 -0.7590619907  0.5219530189
        H 0.0000000000  0.7590619907  0.5219530189
        """)

    @staticmethod
    def create_atomic_input(molecule: Molecule, method: str = "hf", basis: str = "sto-3g") -> AtomicInput:
        """Create AtomicInput for QCEngine."""
        return qcel.models.AtomicInput(
            molecule=molecule,
            driver="energy",
            model={"method": method, "basis": basis}
        )

    @staticmethod
    def qcengine_worker_task(
        worker_id: int,
        task_config: Dict,
        molecule: Molecule,
        method: str = "hf",
        basis: str = "sto-3g"
    ) -> Tuple[int, bool, str, float]:
        """
        QCEngine worker task for threading tests.

        Returns:
            Tuple of (worker_id, success, error_message, energy)
        """
        try:
            # Create atomic input
            atomic_input = ThreadingQCEngineValidator.create_atomic_input(molecule, method, basis)

            # Execute with QCEngine
            result = qcng.compute(atomic_input, "psi4", task_config=task_config, raise_error=True)

            if result.success:
                return worker_id, True, "", result.return_result
            else:
                return worker_id, False, str(result.error), 0.0

        except Exception as e:
            return worker_id, False, str(e), 0.0


def test_qcengine_single_thread():
    """Test QCEngine works in single thread (baseline)."""
    validator = ThreadingQCEngineValidator()
    molecule = validator.create_test_molecule()

    task_config = {
        "memory": 1.0,
        "ncores": 1,
        "nnodes": 1,
        "cores_per_rank": 1,
        "retries": 0
    }

    worker_id, success, error, energy = validator.qcengine_worker_task(
        0, task_config, molecule
    )

    assert success, f"Single thread QCEngine failed: {error}"
    assert abs(energy - (-76.026760737428)) < 1e-6, f"Unexpected energy: {energy}"


def test_qcengine_threading_basic():
    """Test QCEngine with basic threading (2 workers)."""
    validator = ThreadingQCEngineValidator()
    molecule = validator.create_test_molecule()

    task_config = {
        "memory": 1.0,
        "ncores": 1,
        "nnodes": 1,
        "cores_per_rank": 1,
        "retries": 0
    }

    results = []
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(validator.qcengine_worker_task, i, task_config, molecule)
            for i in range(2)
        ]

        for future in as_completed(futures):
            results.append(future.result())

    # Validate all workers succeeded
    for worker_id, success, error, energy in results:
        assert success, f"Worker {worker_id} failed: {error}"
        assert abs(energy - (-76.026760737428)) < 1e-6, f"Worker {worker_id} unexpected energy: {energy}"


def test_qcengine_threading_stress():
    """Test QCEngine with stress threading (4 concurrent workers)."""
    validator = ThreadingQCEngineValidator()
    molecule = validator.create_test_molecule()

    task_config = {
        "memory": 1.0,
        "ncores": 1,
        "nnodes": 1,
        "cores_per_rank": 1,
        "retries": 0
    }

    num_workers = 4
    results = []

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(validator.qcengine_worker_task, i, task_config, molecule)
            for i in range(num_workers)
        ]

        for future in as_completed(futures):
            results.append(future.result())

    # Validate all workers succeeded
    assert len(results) == num_workers, f"Expected {num_workers} results, got {len(results)}"

    for worker_id, success, error, energy in results:
        assert success, f"Worker {worker_id} failed: {error}"
        assert abs(energy - (-76.026760737428)) < 1e-6, f"Worker {worker_id} unexpected energy: {energy}"


def test_qcengine_threading_different_molecules():
    """Test QCEngine threading with different molecules (more realistic scenario)."""
    validator = ThreadingQCEngineValidator()

    # Create different test molecules
    molecules = [
        qcel.models.Molecule.from_data("H 0 0 0\nH 0 0 1.5"),  # H2
        qcel.models.Molecule.from_data("He 0 0 0"),  # He
        validator.create_test_molecule(),  # H2O
    ]

    expected_energies = [
        -1.1170045598,  # H2/HF/STO-3G
        -2.8351061495,  # He/HF/STO-3G
        -76.026760737428,  # H2O/HF/STO-3G
    ]

    task_config = {
        "memory": 1.0,
        "ncores": 1,
        "nnodes": 1,
        "cores_per_rank": 1,
        "retries": 0
    }

    results = []
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [
            executor.submit(validator.qcengine_worker_task, i, task_config, molecules[i])
            for i in range(len(molecules))
        ]

        for future in as_completed(futures):
            results.append(future.result())

    # Validate all workers succeeded with correct energies
    results.sort(key=lambda x: x[0])  # Sort by worker_id

    for i, (worker_id, success, error, energy) in enumerate(results):
        assert success, f"Worker {worker_id} failed: {error}"
        assert abs(energy - expected_energies[i]) < 1e-6, \
            f"Worker {worker_id} unexpected energy: {energy}, expected: {expected_energies[i]}"


def test_qcengine_threading_config_isolation():
    """Test that different task_config values work correctly in parallel."""
    validator = ThreadingQCEngineValidator()
    molecule = validator.create_test_molecule()

    # Different memory configurations
    task_configs = [
        {"memory": 0.5, "ncores": 1, "nnodes": 1, "cores_per_rank": 1, "retries": 0},
        {"memory": 1.0, "ncores": 1, "nnodes": 1, "cores_per_rank": 1, "retries": 0},
        {"memory": 1.5, "ncores": 1, "nnodes": 1, "cores_per_rank": 1, "retries": 0},
    ]

    results = []
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [
            executor.submit(validator.qcengine_worker_task, i, task_configs[i], molecule)
            for i in range(len(task_configs))
        ]

        for future in as_completed(futures):
            results.append(future.result())

    # All should succeed despite different configurations
    for worker_id, success, error, energy in results:
        assert success, f"Worker {worker_id} with config {task_configs[worker_id]} failed: {error}"
        assert abs(energy - (-76.026760737428)) < 1e-6, f"Worker {worker_id} unexpected energy: {energy}"


@pytest.mark.slow
def test_qcengine_threading_repeated_execution():
    """Test repeated threading execution to check for resource leaks."""
    validator = ThreadingQCEngineValidator()
    molecule = validator.create_test_molecule()

    task_config = {
        "memory": 1.0,
        "ncores": 1,
        "nnodes": 1,
        "cores_per_rank": 1,
        "retries": 0
    }

    # Run multiple rounds of threading
    for round_num in range(3):
        results = []
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [
                executor.submit(validator.qcengine_worker_task, i, task_config, molecule)
                for i in range(2)
            ]

            for future in as_completed(futures):
                results.append(future.result())

        # Validate all workers succeeded in this round
        for worker_id, success, error, energy in results:
            assert success, f"Round {round_num}, Worker {worker_id} failed: {error}"
            assert abs(energy - (-76.026760737428)) < 1e-6, \
                f"Round {round_num}, Worker {worker_id} unexpected energy: {energy}"


if __name__ == "__main__":
    # Run basic tests if executed directly
    if HAS_QCENGINE and HAS_PSI4:
        print("Running QCEngine threading tests...")

        try:
            test_qcengine_single_thread()
            print("✓ Single thread test passed")
        except Exception as e:
            print(f"✗ Single thread test failed: {e}")

        try:
            test_qcengine_threading_basic()
            print("✓ Basic threading test passed")
        except Exception as e:
            print(f"✗ Basic threading test failed: {e}")

        try:
            test_qcengine_threading_stress()
            print("✓ Stress threading test passed")
        except Exception as e:
            print(f"✗ Stress threading test failed: {e}")

    else:
        print("QCEngine or Psi4 not available - skipping tests")