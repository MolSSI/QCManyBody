"""
Checkpointing infrastructure for parallel execution.

This module provides functionality to save and restore computation state,
enabling resumption of interrupted calculations and incremental result caching.
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Set
from dataclasses import dataclass, asdict
from datetime import datetime

from .task import ParallelTask, TaskResult, TaskStatus

logger = logging.getLogger(__name__)


@dataclass
class CheckpointMetadata:
    """
    Metadata for checkpoint file.

    Attributes
    ----------
    created_at : str
        ISO timestamp of checkpoint creation
    qcmanybody_version : str
        QCManyBody version
    total_tasks : int
        Total number of tasks in computation
    completed_tasks : int
        Number of completed tasks
    failed_tasks : int
        Number of failed tasks
    in_progress_tasks : int
        Number of tasks in progress (not yet finished)
    """

    created_at: str
    qcmanybody_version: str
    total_tasks: int
    completed_tasks: int
    failed_tasks: int
    in_progress_tasks: int


class CheckpointManager:
    """
    Manages checkpointing for parallel execution.

    The CheckpointManager handles saving and loading of computation state,
    enabling:
    - Resumption of interrupted calculations
    - Incremental result saving
    - Fault tolerance

    **Checkpoint File Format:**

    Checkpoint files are JSON with the following structure:
    ```json
    {
      "metadata": {
        "created_at": "2025-01-15T10:30:00",
        "qcmanybody_version": "0.1.0",
        "total_tasks": 100,
        "completed_tasks": 45,
        "failed_tasks": 2,
        "in_progress_tasks": 0
      },
      "results": {
        "task_id_1": {...},
        "task_id_2": {...}
      }
    }
    ```

    Parameters
    ----------
    checkpoint_file : str or Path
        Path to checkpoint file
    auto_save : bool
        Automatically save after each result (default: True)
    save_interval : int
        Save every N results (if auto_save=True)

    Examples
    --------
    >>> # Create checkpoint manager
    >>> manager = CheckpointManager("checkpoint.json", auto_save=True, save_interval=10)
    >>>
    >>> # Check if checkpoint exists
    >>> if manager.exists():
    ...     # Load existing results
    ...     results = manager.load()
    ...     print(f"Resuming from {len(results)} completed tasks")
    >>>
    >>> # Save result
    >>> manager.save_result(result)
    >>>
    >>> # Get completed task IDs
    >>> completed_ids = manager.get_completed_task_ids()
    >>>
    >>> # Finalize checkpoint
    >>> all_results = manager.finalize()
    """

    def __init__(
        self,
        checkpoint_file: str | Path,
        auto_save: bool = True,
        save_interval: int = 10
    ):
        """
        Initialize checkpoint manager.

        Parameters
        ----------
        checkpoint_file : str or Path
            Path to checkpoint file
        auto_save : bool
            Automatically save after each result
        save_interval : int
            Save every N results (if auto_save=True)
        """
        self.checkpoint_file = Path(checkpoint_file)
        self.auto_save = auto_save
        self.save_interval = save_interval

        self._results: Dict[str, TaskResult] = {}
        self._save_counter = 0
        self._metadata: Optional[CheckpointMetadata] = None

    def exists(self) -> bool:
        """
        Check if checkpoint file exists.

        Returns
        -------
        bool
            True if checkpoint file exists
        """
        return self.checkpoint_file.exists()

    def load(self) -> List[TaskResult]:
        """
        Load results from checkpoint file.

        Returns
        -------
        List[TaskResult]
            Previously computed results

        Raises
        ------
        FileNotFoundError
            If checkpoint file doesn't exist
        ValueError
            If checkpoint file is invalid
        """
        if not self.exists():
            raise FileNotFoundError(f"Checkpoint file not found: {self.checkpoint_file}")

        logger.info(f"Loading checkpoint from {self.checkpoint_file}")

        try:
            with open(self.checkpoint_file, 'r') as f:
                data = json.load(f)

            # Load metadata
            if "metadata" in data:
                self._metadata = CheckpointMetadata(**data["metadata"])
                logger.info(
                    f"Checkpoint metadata: {self._metadata.completed_tasks} completed, "
                    f"{self._metadata.failed_tasks} failed out of "
                    f"{self._metadata.total_tasks} total tasks"
                )

            # Load results
            results = []
            if "results" in data:
                for task_id, result_dict in data["results"].items():
                    result = TaskResult(
                        task_id=result_dict["task_id"],
                        success=result_dict["success"],
                        error_message=result_dict.get("error_message"),
                        atomic_result=result_dict.get("atomic_result"),
                        execution_time=result_dict.get("execution_time", 0.0),
                        worker_id=result_dict.get("worker_id")
                    )
                    results.append(result)
                    self._results[task_id] = result

            logger.info(f"Loaded {len(results)} results from checkpoint")
            return results

        except Exception as e:
            raise ValueError(f"Failed to load checkpoint: {str(e)}") from e

    def save_result(self, result: TaskResult) -> None:
        """
        Save a single result to checkpoint.

        Parameters
        ----------
        result : TaskResult
            Result to save
        """
        self._results[result.task_id] = result
        self._save_counter += 1

        if self.auto_save and self._save_counter >= self.save_interval:
            self.save()
            self._save_counter = 0

    def save_results(self, results: List[TaskResult]) -> None:
        """
        Save multiple results to checkpoint.

        Parameters
        ----------
        results : List[TaskResult]
            Results to save
        """
        for result in results:
            self._results[result.task_id] = result

        if self.auto_save:
            self.save()
            self._save_counter = 0

    def save(self, total_tasks: Optional[int] = None) -> None:
        """
        Save current state to checkpoint file.

        Parameters
        ----------
        total_tasks : int, optional
            Total number of tasks (if known)
        """
        # Count task statuses
        completed = sum(1 for r in self._results.values() if r.success)
        failed = sum(1 for r in self._results.values() if not r.success)
        total = total_tasks or len(self._results)

        # Create metadata
        metadata = CheckpointMetadata(
            created_at=datetime.now().isoformat(),
            qcmanybody_version="0.1.0",  # TODO: Get from package
            total_tasks=total,
            completed_tasks=completed,
            failed_tasks=failed,
            in_progress_tasks=max(0, total - completed - failed)
        )

        # Prepare checkpoint data
        checkpoint_data = {
            "metadata": asdict(metadata),
            "results": {
                task_id: self._result_to_dict(result)
                for task_id, result in self._results.items()
            }
        }

        # Ensure directory exists
        self.checkpoint_file.parent.mkdir(parents=True, exist_ok=True)

        # Write checkpoint (atomic write using temp file)
        temp_file = self.checkpoint_file.with_suffix('.tmp')
        try:
            with open(temp_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)

            # Atomic rename
            temp_file.replace(self.checkpoint_file)

            logger.debug(
                f"Saved checkpoint: {completed} completed, {failed} failed "
                f"out of {total} total tasks"
            )

        except Exception as e:
            logger.error(f"Failed to save checkpoint: {str(e)}")
            if temp_file.exists():
                temp_file.unlink()
            raise

    def get_completed_task_ids(self) -> Set[str]:
        """
        Get IDs of completed tasks.

        Returns
        -------
        Set[str]
            Set of completed task IDs
        """
        return {
            task_id for task_id, result in self._results.items()
            if result.success
        }

    def get_failed_task_ids(self) -> Set[str]:
        """
        Get IDs of failed tasks.

        Returns
        -------
        Set[str]
            Set of failed task IDs
        """
        return {
            task_id for task_id, result in self._results.items()
            if not result.success
        }

    def get_result(self, task_id: str) -> Optional[TaskResult]:
        """
        Get result for specific task.

        Parameters
        ----------
        task_id : str
            Task ID

        Returns
        -------
        TaskResult or None
            Result if found, None otherwise
        """
        return self._results.get(task_id)

    def has_result(self, task_id: str) -> bool:
        """
        Check if result exists for task.

        Parameters
        ----------
        task_id : str
            Task ID

        Returns
        -------
        bool
            True if result exists
        """
        return task_id in self._results

    def filter_pending_tasks(self, tasks: List[ParallelTask]) -> List[ParallelTask]:
        """
        Filter out tasks that already have results.

        This is useful for resuming calculations - only run tasks
        that haven't been completed yet.

        Parameters
        ----------
        tasks : List[ParallelTask]
            All tasks

        Returns
        -------
        List[ParallelTask]
            Tasks without results (pending tasks)

        Examples
        --------
        >>> manager = CheckpointManager("checkpoint.json")
        >>> if manager.exists():
        ...     manager.load()
        ...     tasks = manager.filter_pending_tasks(all_tasks)
        ...     print(f"Resuming: {len(tasks)} tasks remaining")
        """
        completed_ids = self.get_completed_task_ids()
        pending = [task for task in tasks if task.task_id not in completed_ids]

        if pending:
            logger.info(
                f"Filtered tasks: {len(completed_ids)} completed, "
                f"{len(pending)} pending"
            )

        return pending

    def finalize(self) -> List[TaskResult]:
        """
        Finalize checkpoint and return all results.

        Performs final save and returns all accumulated results.

        Returns
        -------
        List[TaskResult]
            All results
        """
        # Final save
        if self._results:
            self.save()

        logger.info(f"Finalized checkpoint with {len(self._results)} results")
        return list(self._results.values())

    def clear(self) -> None:
        """
        Clear checkpoint data.

        Removes in-memory results but does not delete checkpoint file.
        """
        self._results.clear()
        self._save_counter = 0
        logger.debug("Cleared checkpoint data")

    def delete(self) -> None:
        """
        Delete checkpoint file.

        Removes checkpoint file from disk and clears in-memory data.
        """
        if self.checkpoint_file.exists():
            self.checkpoint_file.unlink()
            logger.info(f"Deleted checkpoint file: {self.checkpoint_file}")

        self.clear()

    @staticmethod
    def _result_to_dict(result: TaskResult) -> Dict[str, Any]:
        """Convert TaskResult to JSON-serializable dict."""
        return {
            "task_id": result.task_id,
            "success": result.success,
            "error_message": result.error_message,
            "atomic_result": result.atomic_result,
            "execution_time": result.execution_time,
            "worker_id": result.worker_id
        }


def create_checkpoint_manager(
    checkpoint_file: str | Path,
    resume_if_exists: bool = True
) -> CheckpointManager:
    """
    Convenience function to create and optionally load checkpoint manager.

    Parameters
    ----------
    checkpoint_file : str or Path
        Path to checkpoint file
    resume_if_exists : bool
        If True and checkpoint exists, load it automatically

    Returns
    -------
    CheckpointManager
        Checkpoint manager (loaded if resume_if_exists=True)

    Examples
    --------
    >>> # Create and auto-load if exists
    >>> manager = create_checkpoint_manager("checkpoint.json", resume_if_exists=True)
    >>> if manager.exists():
    ...     print(f"Resumed from checkpoint")
    """
    manager = CheckpointManager(checkpoint_file)

    if resume_if_exists and manager.exists():
        try:
            manager.load()
            logger.info("Resumed from existing checkpoint")
        except Exception as e:
            logger.warning(f"Failed to load checkpoint: {e}. Starting fresh.")

    return manager
