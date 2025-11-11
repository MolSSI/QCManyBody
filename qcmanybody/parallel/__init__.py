"""
QCManyBody Parallel Execution Module

This module provides parallel execution capabilities for many-body expansion
calculations. It supports multiple execution backends:
- Sequential (reference implementation)
- Multiprocessing (single-node parallelism)
- MPI (distributed multi-node parallelism, optional dependency)

Basic Usage
-----------
>>> from qcmanybody import ManyBodyComputer
>>> result = ManyBodyComputer.from_manybodyinput(
...     input_model,
...     parallel=True,
...     n_workers=4
... )

Advanced Usage
--------------
>>> from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig
>>> config = ExecutorConfig(n_workers=8, timeout_per_task=3600)
>>> executor = MultiprocessingExecutor(config)
>>> result = ManyBodyComputer.from_manybodyinput(
...     input_model,
...     executor=executor
... )
"""

from .base import BaseParallelExecutor, ExecutorConfig
from .task import ParallelTask, TaskResult
from .executors.sequential import SequentialExecutor

__all__ = [
    "BaseParallelExecutor",
    "ExecutorConfig",
    "ParallelTask",
    "TaskResult",
    "SequentialExecutor",
]

# Conditionally import multiprocessing executor
try:
    from .executors.multiprocessing import MultiprocessingExecutor
    __all__.append("MultiprocessingExecutor")
except ImportError:
    pass

# Conditionally import MPI executor
try:
    from .executors.mpi import MPIExecutor
    __all__.append("MPIExecutor")
except (ImportError, ModuleNotFoundError):
    pass

__version__ = "0.1.0"
