"""
Task and result data models for parallel execution.

This module defines the dataclasses used to specify computational tasks
and store their results.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List
from enum import Enum
import time


class TaskStatus(Enum):
    """Status of a parallel task."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    TIMEOUT = "timeout"
    CANCELLED = "cancelled"


@dataclass
class ParallelTask:
    """Specification for a single parallel quantum chemistry computation.

    This dataclass encapsulates all information needed to execute a
    single QC calculation in parallel.

    Parameters
    ----------
    task_id : str
        Unique identifier for this task (typically from labeler)
    chemistry : str
        Model chemistry string (e.g., "mp2/cc-pvdz")
    label : str
        Full label string from labeler
    molecule : Any
        QCElemental Molecule object (type hint Any to avoid import)
    atomic_input : Any
        QCElemental AtomicInput object
    priority : int
        Task priority (higher = execute sooner). Default: 0
    estimated_cost : float
        Relative computational cost estimate. Default: 1.0
    nbody : int
        N-body level of this calculation. Default: 1
    depends_on : List[str]
        Task IDs that must complete before this task. Default: []
    metadata : Dict[str, Any]
        Additional task metadata. Default: {}

    Examples
    --------
    >>> from qcelemental.models import Molecule, AtomicInput
    >>> mol = Molecule(symbols=["He", "He"], geometry=[...])
    >>> inp = AtomicInput(molecule=mol, model={...}, driver="energy")
    >>> task = ParallelTask(
    ...     task_id="he2_dimer",
    ...     chemistry="mp2/cc-pvdz",
    ...     label='["mp2/cc-pvdz", [1, 2], [1, 2]]',
    ...     molecule=mol,
    ...     atomic_input=inp,
    ...     priority=1,
    ...     nbody=2
    ... )
    """

    task_id: str
    chemistry: str
    label: str
    molecule: Any  # QCElemental Molecule
    atomic_input: Any  # QCElemental AtomicInput

    # Scheduling metadata
    priority: int = 0
    estimated_cost: float = 1.0
    nbody: int = 1
    depends_on: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __hash__(self):
        """Hash based on task_id for use in sets/dicts."""
        return hash(self.task_id)

    def __eq__(self, other):
        """Equality based on task_id."""
        if not isinstance(other, ParallelTask):
            return False
        return self.task_id == other.task_id

    def __lt__(self, other):
        """Comparison for priority queue (higher priority first)."""
        if not isinstance(other, ParallelTask):
            return NotImplemented
        # Higher priority comes first, then by estimated cost (larger first)
        return (self.priority, self.estimated_cost) > (other.priority, other.estimated_cost)

    def __repr__(self):
        """String representation for debugging."""
        return (
            f"ParallelTask(task_id='{self.task_id}', "
            f"chemistry='{self.chemistry}', nbody={self.nbody}, "
            f"priority={self.priority}, cost={self.estimated_cost:.1f})"
        )


@dataclass
class TaskResult:
    """Result of a parallel task execution.

    This dataclass stores the outcome of executing a ParallelTask,
    including success/failure status, QC results, errors, and timing.

    Parameters
    ----------
    task_id : str
        ID of the task that produced this result
    success : bool
        True if task completed successfully
    status : TaskStatus
        Final status of the task. Default: COMPLETED
    atomic_result : Optional[Any]
        QCElemental AtomicResult object (if successful). Default: None
    error_type : Optional[str]
        Exception type if failed. Default: None
    error_message : Optional[str]
        Error message if failed. Default: None
    error_traceback : Optional[str]
        Full traceback if failed. Default: None
    execution_time : float
        Task execution time in seconds. Default: 0.0
    queue_time : float
        Time spent waiting in queue in seconds. Default: 0.0
    attempt_number : int
        Attempt number (for retries). Default: 1
    worker_id : Optional[str]
        ID of worker that executed this task. Default: None
    timestamp : float
        Unix timestamp when result was created. Default: current time
    metadata : Dict[str, Any]
        Additional result metadata. Default: {}

    Examples
    --------
    >>> result = TaskResult(
    ...     task_id="he2_dimer",
    ...     success=True,
    ...     atomic_result=atomic_result_obj,
    ...     execution_time=45.2,
    ...     worker_id="worker_0"
    ... )
    >>> if result.success:
    ...     energy = result.return_result
    """

    task_id: str
    success: bool
    status: TaskStatus = TaskStatus.COMPLETED

    # QC result (if successful)
    atomic_result: Optional[Any] = None  # QCElemental AtomicResult

    # Error info (if failed)
    error_type: Optional[str] = None
    error_message: Optional[str] = None
    error_traceback: Optional[str] = None

    # Timing metadata
    execution_time: float = 0.0
    queue_time: float = 0.0
    attempt_number: int = 1
    worker_id: Optional[str] = None
    timestamp: float = field(default_factory=time.time)

    # Additional metadata
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def return_result(self) -> Optional[Any]:
        """Extract primary result (energy/gradient/hessian) from atomic_result.

        Returns
        -------
        Optional[Any]
            The return_result from AtomicResult, or None if not available
        """
        if self.success and self.atomic_result is not None:
            return getattr(self.atomic_result, "return_result", None)
        return None

    @property
    def total_time(self) -> float:
        """Total time including queue and execution.

        Returns
        -------
        float
            queue_time + execution_time
        """
        return self.queue_time + self.execution_time

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary (for serialization).

        Returns
        -------
        Dict[str, Any]
            Dictionary representation of this result
        """
        return {
            "task_id": self.task_id,
            "success": self.success,
            "status": self.status.value,
            "error_type": self.error_type,
            "error_message": self.error_message,
            "execution_time": self.execution_time,
            "queue_time": self.queue_time,
            "attempt_number": self.attempt_number,
            "worker_id": self.worker_id,
            "timestamp": self.timestamp,
            "metadata": self.metadata,
            # Note: atomic_result is complex and not included in simple dict
        }

    def __repr__(self):
        """String representation for debugging."""
        status_symbol = "✓" if self.success else "✗"
        return (
            f"TaskResult({status_symbol} task_id='{self.task_id}', "
            f"status={self.status.value}, time={self.execution_time:.2f}s, "
            f"attempt={self.attempt_number})"
        )
