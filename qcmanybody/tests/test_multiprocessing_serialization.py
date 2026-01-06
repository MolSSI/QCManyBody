"""Regression tests for multiprocessing execution under spawn start method."""

import importlib.util
import multiprocessing as mp
from pathlib import Path

import pytest

from qcmanybody.parallel import ParallelConfig, ParallelManyBodyExecutor

_SCRIPT_PATH = (
    Path(__file__).resolve().parents[2]
    / "parallel-execution-project"
    / "tests"
    / "test_water4_mbe4_multiprocessing.py"
)


def _load_water4_module():
    spec = importlib.util.spec_from_file_location("water4_multiprocessing_module", _SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load water4 multiprocessing module from {_SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def water4_module():
    return _load_water4_module()


@pytest.fixture(scope="module")
def ensure_spawn_start_method():
    start_method = mp.get_start_method(allow_none=True)
    if start_method != "spawn":
        try:
            mp.set_start_method("spawn", force=True)
        except RuntimeError as exc:  # pragma: no cover - system-dependent guard
            pytest.skip(f"Unable to force spawn start method: {exc}")
    return "spawn"


def test_multiprocessing_spawn_execution_matches_serial(water4_module, ensure_spawn_start_method):
    molecule_serial = water4_module.build_water4_molecule()
    molecule_parallel = water4_module.build_water4_molecule()

    serial_core = water4_module.build_manybody_core(molecule_serial)
    parallel_core = water4_module.build_manybody_core(molecule_parallel)

    serial_config = ParallelConfig(max_workers=1, execution_mode="serial", use_qcengine=False)
    serial_executor = ParallelManyBodyExecutor(serial_core, serial_config)
    serial_results = serial_executor.execute_full_calculation()

    multiprocessing_config = ParallelConfig(max_workers=2, execution_mode="multiprocessing", use_qcengine=False)
    multiprocessing_executor = ParallelManyBodyExecutor(parallel_core, multiprocessing_config)
    multiprocessing_results = multiprocessing_executor.execute_full_calculation()

    serial_executor.validate_parallel_correctness(multiprocessing_results, serial_results)

    assert multiprocessing_results
    assert set(multiprocessing_results.keys()) == set(serial_results.keys())