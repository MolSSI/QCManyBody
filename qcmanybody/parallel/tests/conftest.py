"""
Pytest fixtures for parallel module tests.
"""

import pytest
import numpy as np


@pytest.fixture
def helium_dimer():
    """Create a simple helium dimer molecule."""
    from qcelemental.models import Molecule

    return Molecule(
        symbols=["He", "He"],
        geometry=[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 3.0],
        ],
        fragments=[[0], [1]],
    )


@pytest.fixture
def helium_trimer():
    """Create a helium trimer molecule."""
    from qcelemental.models import Molecule

    return Molecule(
        symbols=["He", "He", "He"],
        geometry=[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 3.0],
            [0.0, 3.0, 0.0],
        ],
        fragments=[[0], [1], [2]],
    )


@pytest.fixture
def water_dimer():
    """Create a water dimer molecule."""
    from qcelemental.models import Molecule

    return Molecule(
        **{
            "symbols": ["O", "H", "H", "O", "H", "H"],
            "geometry": [
                [0.0, 0.0, 0.0],
                [0.0, 1.5, 0.0],
                [1.5, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [5.0, 1.5, 0.0],
                [6.5, 0.0, 0.0],
            ],
            "fragments": [[0, 1, 2], [3, 4, 5]],
        }
    )


@pytest.fixture
def mock_atomic_input(helium_dimer):
    """Create a mock AtomicInput for testing."""
    from qcelemental.models import AtomicInput

    return AtomicInput(
        molecule=helium_dimer,
        driver="energy",
        model={"method": "hf", "basis": "sto-3g"},
        keywords={},
    )


@pytest.fixture
def mock_parallel_task(mock_atomic_input):
    """Create a mock ParallelTask for testing."""
    from qcmanybody.parallel import ParallelTask

    return ParallelTask(
        task_id="test_task_1",
        chemistry="hf/sto-3g",
        label='["hf/sto-3g", [1], [1]]',
        molecule=mock_atomic_input.molecule,
        atomic_input=mock_atomic_input,
        priority=0,
        nbody=1,
    )


@pytest.fixture
def executor_config():
    """Create a default ExecutorConfig for testing."""
    from qcmanybody.parallel import ExecutorConfig

    return ExecutorConfig(
        n_workers=2,
        timeout_per_task=60.0,  # Short timeout for tests
        max_retries=1,
        checkpoint_interval=5,
        log_level="DEBUG",
    )


@pytest.fixture
def sequential_executor(executor_config):
    """Create a SequentialExecutor for testing."""
    from qcmanybody.parallel import SequentialExecutor

    return SequentialExecutor(executor_config)


@pytest.fixture
def multiprocessing_executor(executor_config):
    """Create a MultiprocessingExecutor for testing."""
    pytest.importorskip("qcengine")
    from qcmanybody.parallel.executors.multiprocessing import MultiprocessingExecutor

    return MultiprocessingExecutor(executor_config)
