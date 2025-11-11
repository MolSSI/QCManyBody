"""
Base classes and interfaces for parallel execution.

This module defines the abstract interface that all parallel executors
must implement, along with configuration dataclasses.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Optional, Callable, Any, Dict
import logging

logger = logging.getLogger(__name__)


@dataclass
class ExecutorConfig:
    """Configuration for parallel executors.

    Parameters
    ----------
    n_workers : Optional[int]
        Number of parallel workers. If None, auto-detect based on
        available CPU cores. Default: None
    timeout_per_task : float
        Maximum time (in seconds) allowed for each task. Default: 3600.0
    max_retries : int
        Maximum number of retry attempts for failed tasks. Default: 2
    checkpoint_interval : int
        Save checkpoint after this many completed tasks. Default: 10
    checkpoint_file : Optional[str]
        Path to checkpoint file. If None, no checkpointing. Default: None
    cache_dir : Optional[str]
        Directory for result caching. If None, no caching. Default: None
    log_level : str
        Logging level for executor. Default: "INFO"
    scratch_dir : Optional[str]
        Scratch directory for temporary files. Default: None
    """

    n_workers: Optional[int] = None
    timeout_per_task: float = 3600.0
    max_retries: int = 2
    checkpoint_interval: int = 10
    checkpoint_file: Optional[str] = None
    cache_dir: Optional[str] = None
    log_level: str = "INFO"
    scratch_dir: Optional[str] = None

    def __post_init__(self):
        """Validate configuration parameters."""
        if self.n_workers is not None and self.n_workers < 1:
            raise ValueError(f"n_workers must be >= 1, got {self.n_workers}")
        if self.timeout_per_task <= 0:
            raise ValueError(f"timeout_per_task must be > 0, got {self.timeout_per_task}")
        if self.max_retries < 0:
            raise ValueError(f"max_retries must be >= 0, got {self.max_retries}")
        if self.checkpoint_interval < 1:
            raise ValueError(f"checkpoint_interval must be >= 1, got {self.checkpoint_interval}")


class BaseParallelExecutor(ABC):
    """Abstract base class for parallel executors.

    All executors must implement this interface to be compatible with
    ParallelManyBodyComputer. Executors are responsible for:
    1. Managing worker processes/threads/nodes
    2. Distributing tasks to workers
    3. Collecting and returning results
    4. Handling errors and retries
    5. Resource cleanup

    Executors can be used as context managers::

        with MyExecutor(config) as executor:
            results = executor.execute(tasks)

    Parameters
    ----------
    config : Optional[ExecutorConfig]
        Executor configuration. If None, uses default configuration.

    Attributes
    ----------
    config : ExecutorConfig
        Configuration for this executor
    """

    def __init__(self, config: Optional[ExecutorConfig] = None):
        self.config = config or ExecutorConfig()
        self._is_initialized = False
        self._setup_logging()

    def _setup_logging(self):
        """Configure logging for this executor."""
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.logger.setLevel(self.config.log_level)

    @abstractmethod
    def initialize(self) -> None:
        """Initialize executor resources.

        This method is called before task execution begins. It should:
        - Create worker pools/processes
        - Initialize communication channels
        - Set up scratch directories
        - Perform any necessary validation

        Raises
        ------
        RuntimeError
            If initialization fails
        """
        pass

    @abstractmethod
    def execute(
        self,
        tasks: List["ParallelTask"],
        progress_callback: Optional[Callable[[str, int, int], None]] = None,
    ) -> List["TaskResult"]:
        """Execute tasks in parallel.

        This is the main execution method. It should:
        1. Distribute tasks to available workers
        2. Monitor task execution
        3. Handle failures and retries
        4. Call progress_callback with updates
        5. Return results in same order as input tasks

        Parameters
        ----------
        tasks : List[ParallelTask]
            List of tasks to execute. Order is preserved in results.
        progress_callback : Optional[Callable[[str, int, int], None]]
            Optional callback function called as:
            progress_callback(task_id, completed_count, total_count)
            Can be used to update progress bars or logging.

        Returns
        -------
        List[TaskResult]
            Results in the same order as input tasks. Failed tasks
            return TaskResult with success=False.

        Raises
        ------
        RuntimeError
            If executor is not initialized
        """
        pass

    @abstractmethod
    def shutdown(self, wait: bool = True) -> None:
        """Shutdown executor and clean up resources.

        This method should:
        - Terminate all worker processes/threads
        - Close communication channels
        - Clean up temporary files
        - Release any held resources

        Parameters
        ----------
        wait : bool
            If True, wait for all workers to finish cleanly.
            If False, terminate workers immediately. Default: True
        """
        pass

    def __enter__(self):
        """Enter context manager - initialize executor."""
        self.initialize()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager - shutdown executor."""
        self.shutdown()
        return False

    @property
    def name(self) -> str:
        """Human-readable name of this executor."""
        return self.__class__.__name__

    @property
    def is_initialized(self) -> bool:
        """Check if executor is initialized and ready."""
        return self._is_initialized

    def validate_tasks(self, tasks: List["ParallelTask"]) -> None:
        """Validate task list before execution.

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to validate

        Raises
        ------
        ValueError
            If task list is invalid
        """
        if not tasks:
            raise ValueError("Task list cannot be empty")

        # Check for duplicate task IDs
        task_ids = [t.task_id for t in tasks]
        if len(task_ids) != len(set(task_ids)):
            duplicates = [tid for tid in task_ids if task_ids.count(tid) > 1]
            raise ValueError(f"Duplicate task IDs found: {set(duplicates)}")

    def get_info(self) -> Dict[str, Any]:
        """Get information about this executor.

        Returns
        -------
        Dict[str, Any]
            Dictionary with executor metadata:
            - name: Executor class name
            - n_workers: Number of workers
            - is_initialized: Initialization status
            - config: Configuration dict
        """
        return {
            "name": self.name,
            "n_workers": self.config.n_workers,
            "is_initialized": self.is_initialized,
            "config": self.config.__dict__,
        }
