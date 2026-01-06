"""
QCManyBody Parallel Execution Module

This module provides parallel execution capabilities for many-body expansion
calculations. It supports multiple execution backends:
- Sequential (reference implementation)
- Multiprocessing (single-node parallelism)
- ConcurrentExecutor (single-node with concurrent.futures)
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
from .task import ParallelTask, TaskResult, TaskStatus
from .executors.sequential import SequentialExecutor
from .computer_parallel import ParallelManyBodyComputer, parallel_compute_from_manybodyinput
from .checkpoint import CheckpointManager, CheckpointMetadata, create_checkpoint_manager
from .scheduler import TaskScheduler, SchedulingStrategy, assign_task_priorities, estimate_load_balance

__all__ = [
    "BaseParallelExecutor",
    "ExecutorConfig",
    "ParallelTask",
    "TaskResult",
    "TaskStatus",
    "SequentialExecutor",
    "ParallelManyBodyComputer",
    "parallel_compute_from_manybodyinput",
    "CheckpointManager",
    "CheckpointMetadata",
    "create_checkpoint_manager",
    "TaskScheduler",
    "SchedulingStrategy",
    "assign_task_priorities",
    "estimate_load_balance",
]

# Conditionally import multiprocessing executor
try:
    from .executors.multiprocessing import MultiprocessingExecutor
    __all__.append("MultiprocessingExecutor")
except ImportError:
    pass

# Conditionally import concurrent.futures executor
try:
    from .executors.concurrent import ConcurrentExecutor, ConcurrentExecutorConfig
    __all__.extend(["ConcurrentExecutor", "ConcurrentExecutorConfig"])
except ImportError:
    pass

# Conditionally import MPI executor
try:
    from .executors.mpi import MPIExecutor
    __all__.append("MPIExecutor")
except (ImportError, ModuleNotFoundError):
    pass

__version__ = "0.1.0"
