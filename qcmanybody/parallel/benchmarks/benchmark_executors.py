"""
Performance benchmark for parallel executors.

Compares execution time between Sequential and Multiprocessing executors.
"""

import time
import sys
from typing import List
from unittest.mock import Mock

# Add parent directory to path if running as script
if __name__ == '__main__':
    sys.path.insert(0, '/home/user/QCManyBody')

from qcmanybody.parallel import (
    SequentialExecutor,
    ExecutorConfig,
    ParallelTask,
)
from qcmanybody.parallel.executors.multiprocessing import MultiprocessingExecutor
from qcelemental.models import Molecule, AtomicInput


def create_mock_tasks(n_tasks: int, delay: float = 0.1) -> List[ParallelTask]:
    """Create mock tasks for benchmarking."""
    mol = Molecule(symbols=["He"], geometry=[[0, 0, 0]])
    inp = AtomicInput(molecule=mol, driver="energy", model={"method": "hf", "basis": "sto-3g"})

    tasks = []
    for i in range(n_tasks):
        task = ParallelTask(
            task_id=f"task_{i}",
            chemistry="hf/sto-3g",
            label=f"task_{i}",
            molecule=mol,
            atomic_input=inp,
        )
        tasks.append(task)

    return tasks


def benchmark_sequential(tasks: List[ParallelTask], config: ExecutorConfig) -> float:
    """Benchmark sequential executor."""
    from unittest.mock import patch

    # Mock qcengine.compute to simulate work
    def mock_compute(*args, **kwargs):
        time.sleep(0.1)  # Simulate 100ms computation
        result = Mock()
        result.success = True
        result.return_result = -5.0
        return result

    with patch('qcengine.compute', side_effect=mock_compute):
        executor = SequentialExecutor(config)
        with executor:
            start = time.time()
            results = executor.execute(tasks)
            elapsed = time.time() - start

    return elapsed


def benchmark_multiprocessing(tasks: List[ParallelTask], config: ExecutorConfig) -> float:
    """Benchmark multiprocessing executor."""
    from unittest.mock import patch

    # Mock qcengine.compute to simulate work
    def mock_compute(*args, **kwargs):
        time.sleep(0.1)  # Simulate 100ms computation
        result = Mock()
        result.success = True
        result.return_result = -5.0
        return result

    with patch('qcengine.compute', side_effect=mock_compute):
        executor = MultiprocessingExecutor(config)
        with executor:
            start = time.time()
            results = executor.execute(tasks)
            elapsed = time.time() - start

    return elapsed


def run_benchmarks():
    """Run benchmarks and print results."""
    print("QCManyBody Parallel Execution Benchmarks")
    print("=" * 60)
    print()

    # Test configurations
    test_cases = [
        (4, 2),   # 4 tasks, 2 workers
        (8, 2),   # 8 tasks, 2 workers
        (8, 4),   # 8 tasks, 4 workers
        (16, 4),  # 16 tasks, 4 workers
    ]

    for n_tasks, n_workers in test_cases:
        print(f"Test: {n_tasks} tasks, {n_workers} workers")
        print("-" * 60)

        tasks = create_mock_tasks(n_tasks)
        config = ExecutorConfig(n_workers=n_workers, timeout_per_task=60.0)

        # Benchmark sequential
        seq_time = benchmark_sequential(tasks, ExecutorConfig(n_workers=1))
        print(f"  Sequential:      {seq_time:.3f}s")

        # Benchmark multiprocessing
        mp_time = benchmark_multiprocessing(tasks, config)
        print(f"  Multiprocessing: {mp_time:.3f}s")

        # Calculate speedup
        speedup = seq_time / mp_time if mp_time > 0 else 0
        efficiency = (speedup / n_workers) * 100 if n_workers > 0 else 0

        print(f"  Speedup:         {speedup:.2f}x")
        print(f"  Efficiency:      {efficiency:.1f}%")
        print(f"  Expected:        ~{min(n_workers, n_tasks)}x (theoretical max)")
        print()

    print("=" * 60)
    print("Notes:")
    print("- Each task simulates 100ms of computation")
    print("- Speedup = Sequential Time / Parallel Time")
    print("- Efficiency = (Speedup / Workers) * 100%")
    print("- Actual speedup depends on system cores and overhead")
    print()


if __name__ == '__main__':
    # Check dependencies
    try:
        import qcengine
    except ImportError:
        print("WARNING: qcengine not installed, using mocked compute")
        print()

    run_benchmarks()
