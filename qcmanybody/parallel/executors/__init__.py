"""
Parallel executor implementations.

This package contains concrete implementations of BaseParallelExecutor:
- SequentialExecutor: Reference implementation (no parallelism)
- MultiprocessingExecutor: Single-node parallelism using multiprocessing
- ConcurrentExecutor: Single-node parallelism using concurrent.futures
- MPIExecutor: Multi-node distributed parallelism using MPI (optional)
"""

from .sequential import SequentialExecutor

__all__ = ["SequentialExecutor"]

# Conditionally import multiprocessing executor
try:
    from .multiprocessing import MultiprocessingExecutor

    __all__.append("MultiprocessingExecutor")
except ImportError:
    pass

# Conditionally import concurrent.futures executor
try:
    from .concurrent import ConcurrentExecutor, ConcurrentExecutorConfig

    __all__.extend(["ConcurrentExecutor", "ConcurrentExecutorConfig"])
except ImportError:
    pass

# Conditionally import MPI executor
try:
    from .mpi import MPIExecutor

    __all__.append("MPIExecutor")
except (ImportError, ModuleNotFoundError):
    pass
