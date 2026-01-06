"""
Sequential executor - reference implementation with no parallelism.

This executor runs tasks one at a time in the current process/thread.
It serves as:
1. Reference implementation for correctness testing
2. Fallback when parallelism is unavailable or undesired
3. Debugging tool for isolating parallel vs. serial issues
"""

import logging
from typing import List, Optional, Callable

from ..base import BaseParallelExecutor
from ..task import ParallelTask, TaskResult, TaskStatus
from ..worker import execute_single_task

logger = logging.getLogger(__name__)


class SequentialExecutor(BaseParallelExecutor):
    """Sequential executor - executes tasks one at a time.

    This executor provides no parallelism but ensures correct
    serial execution. Useful for debugging, testing, and as a
    fallback option.

    Examples
    --------
    >>> from qcmanybody.parallel import SequentialExecutor, ExecutorConfig
    >>> config = ExecutorConfig(timeout_per_task=1800)
    >>> executor = SequentialExecutor(config)
    >>> with executor:
    ...     results = executor.execute(tasks)

    Notes
    -----
    This executor ignores the n_workers configuration parameter
    since it always uses exactly 1 worker (the current process).
    """

    def initialize(self) -> None:
        """Initialize sequential executor.

        No actual initialization needed, but we validate configuration
        and set initialized flag.
        """
        self.logger.info("Initializing SequentialExecutor")

        if self.config.n_workers and self.config.n_workers != 1:
            self.logger.warning(
                f"SequentialExecutor ignores n_workers={self.config.n_workers}, "
                "always uses 1 worker"
            )

        self._is_initialized = True
        self.logger.info("SequentialExecutor initialized successfully")

    def execute(
        self,
        tasks: List[ParallelTask],
        progress_callback: Optional[Callable[[str, int, int], None]] = None,
    ) -> List[TaskResult]:
        """Execute tasks sequentially in order.

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to execute
        progress_callback : Optional[Callable[[str, int, int], None]]
            Progress callback function

        Returns
        -------
        List[TaskResult]
            Results in same order as input tasks

        Raises
        ------
        RuntimeError
            If executor not initialized
        """
        if not self._is_initialized:
            raise RuntimeError("Executor not initialized. Call initialize() first or use context manager.")

        # Validate tasks
        self.validate_tasks(tasks)

        results = []
        total = len(tasks)

        self.logger.info(f"Executing {total} tasks sequentially")

        for i, task in enumerate(tasks):
            self.logger.debug(f"Executing task {i+1}/{total}: {task.task_id}")

            # Execute task
            result = execute_single_task(task, self.config)
            results.append(result)

            # Progress callback
            if progress_callback:
                try:
                    progress_callback(task.task_id, i + 1, total)
                except Exception as e:
                    self.logger.warning(f"Progress callback failed: {e}")

            # Log result
            status_symbol = "✓" if result.success else "✗"
            self.logger.info(
                f"{status_symbol} Task {i+1}/{total} ({task.task_id}): "
                f"{result.execution_time:.2f}s, status={result.status.value}"
            )

            if not result.success:
                self.logger.error(
                    f"Task {task.task_id} failed: {result.error_type}: " f"{result.error_message}"
                )

        # Summary
        successful = sum(1 for r in results if r.success)
        failed = total - successful
        total_time = sum(r.execution_time for r in results)

        self.logger.info(
            f"Sequential execution complete: {successful} succeeded, "
            f"{failed} failed, total time {total_time:.2f}s"
        )

        return results

    def shutdown(self, wait: bool = True) -> None:
        """Shutdown sequential executor.

        No cleanup needed, but we reset the initialized flag.

        Parameters
        ----------
        wait : bool
            Ignored for sequential executor
        """
        self.logger.info("Shutting down SequentialExecutor")
        self._is_initialized = False

    def get_info(self):
        """Get executor information."""
        info = super().get_info()
        info["n_workers"] = 1  # Override to always show 1
        return info
