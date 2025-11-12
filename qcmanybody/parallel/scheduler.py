"""
Task scheduling and load balancing for parallel execution.

This module provides intelligent task scheduling strategies to optimize
parallel execution performance.
"""

from typing import List, Callable, Optional, Dict, Any
from dataclasses import dataclass
import logging

from .task import ParallelTask

logger = logging.getLogger(__name__)


@dataclass
class SchedulingStrategy:
    """Configuration for task scheduling strategy."""

    name: str = "priority_first"
    """Scheduling strategy name: 'priority_first', 'cost_first', 'nbody_first', 'dependency_aware', 'fifo'"""

    enable_load_balancing: bool = True
    """Enable dynamic load balancing"""

    chunk_size: Optional[int] = None
    """Task chunk size for batch submission (None = auto)"""

    reorder_tasks: bool = True
    """Allow reordering tasks for optimization"""


class TaskScheduler:
    """
    Intelligent task scheduler for parallel execution.

    Provides various scheduling strategies to optimize task execution:
    - Priority-based scheduling
    - Cost-aware scheduling
    - N-body level scheduling
    - Dependency-aware scheduling
    - Load balancing

    Examples
    --------
    >>> scheduler = TaskScheduler(strategy="nbody_first")
    >>> scheduled_tasks = scheduler.schedule(tasks)
    >>>
    >>> # With custom strategy
    >>> strategy = SchedulingStrategy(
    ...     name="cost_first",
    ...     enable_load_balancing=True,
    ...     chunk_size=10
    ... )
    >>> scheduler = TaskScheduler(strategy=strategy)
    """

    def __init__(
        self,
        strategy: Optional[SchedulingStrategy] = None,
        n_workers: int = 1
    ):
        """
        Initialize task scheduler.

        Parameters
        ----------
        strategy : SchedulingStrategy, optional
            Scheduling strategy configuration
        n_workers : int
            Number of workers (for load balancing)
        """
        self.strategy = strategy or SchedulingStrategy()
        self.n_workers = n_workers
        self._scheduling_stats: Dict[str, Any] = {}

    def schedule(self, tasks: List[ParallelTask]) -> List[ParallelTask]:
        """
        Schedule tasks according to strategy.

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to schedule

        Returns
        -------
        List[ParallelTask]
            Scheduled tasks (possibly reordered)
        """
        if not tasks:
            return []

        if not self.strategy.reorder_tasks:
            logger.debug("Task reordering disabled, returning original order")
            return tasks

        logger.info(f"Scheduling {len(tasks)} tasks with strategy: {self.strategy.name}")

        # Apply scheduling strategy
        if self.strategy.name == "priority_first":
            scheduled = self._schedule_priority_first(tasks)
        elif self.strategy.name == "cost_first":
            scheduled = self._schedule_cost_first(tasks)
        elif self.strategy.name == "nbody_first":
            scheduled = self._schedule_nbody_first(tasks)
        elif self.strategy.name == "dependency_aware":
            scheduled = self._schedule_dependency_aware(tasks)
        elif self.strategy.name == "fifo":
            scheduled = tasks.copy()
        else:
            logger.warning(f"Unknown strategy '{self.strategy.name}', using FIFO")
            scheduled = tasks.copy()

        # Compute and log statistics
        self._compute_scheduling_stats(tasks, scheduled)

        return scheduled

    def create_chunks(
        self,
        tasks: List[ParallelTask]
    ) -> List[List[ParallelTask]]:
        """
        Split tasks into chunks for batch submission.

        This reduces memory footprint by submitting tasks in batches
        rather than all at once.

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to chunk

        Returns
        -------
        List[List[ParallelTask]]
            List of task chunks
        """
        if not tasks:
            return []

        chunk_size = self._determine_chunk_size(tasks)

        chunks = []
        for i in range(0, len(tasks), chunk_size):
            chunk = tasks[i:i + chunk_size]
            chunks.append(chunk)

        logger.info(f"Split {len(tasks)} tasks into {len(chunks)} chunks of size ~{chunk_size}")

        return chunks

    def get_stats(self) -> Dict[str, Any]:
        """
        Get scheduling statistics.

        Returns
        -------
        Dict[str, Any]
            Scheduling statistics
        """
        return self._scheduling_stats.copy()

    # Private methods

    def _schedule_priority_first(self, tasks: List[ParallelTask]) -> List[ParallelTask]:
        """Schedule by priority (highest first)."""
        return sorted(tasks, key=lambda t: -t.priority)

    def _schedule_cost_first(self, tasks: List[ParallelTask]) -> List[ParallelTask]:
        """Schedule by estimated cost (highest first)."""
        return sorted(tasks, key=lambda t: -t.estimated_cost)

    def _schedule_nbody_first(self, tasks: List[ParallelTask]) -> List[ParallelTask]:
        """Schedule by n-body level (highest first)."""
        return sorted(tasks, key=lambda t: -t.nbody)

    def _schedule_dependency_aware(self, tasks: List[ParallelTask]) -> List[ParallelTask]:
        """
        Schedule respecting dependencies.

        Uses topological sort to ensure dependencies are computed first.
        """
        # Build dependency graph
        task_dict = {t.task_id: t for t in tasks}
        in_degree = {t.task_id: 0 for t in tasks}

        # Count incoming edges (dependencies)
        for task in tasks:
            for dep_id in task.depends_on:
                if dep_id in in_degree:
                    in_degree[task.task_id] += 1

        # Topological sort using Kahn's algorithm
        scheduled = []
        queue = [t for t in tasks if in_degree[t.task_id] == 0]

        while queue:
            # Sort queue by priority for tie-breaking
            queue.sort(key=lambda t: -t.priority)

            task = queue.pop(0)
            scheduled.append(task)

            # Update in-degrees of dependent tasks
            for other_task in tasks:
                if task.task_id in other_task.depends_on:
                    in_degree[other_task.task_id] -= 1
                    if in_degree[other_task.task_id] == 0:
                        queue.append(other_task)

        # Check for cycles
        if len(scheduled) != len(tasks):
            logger.warning("Dependency cycle detected, falling back to priority scheduling")
            return self._schedule_priority_first(tasks)

        return scheduled

    def _determine_chunk_size(self, tasks: List[ParallelTask]) -> int:
        """
        Determine optimal chunk size.

        Factors considered:
        - Number of workers
        - Total number of tasks
        - User-specified chunk size
        """
        if self.strategy.chunk_size is not None:
            return self.strategy.chunk_size

        # Auto-determine based on workers and tasks
        n_tasks = len(tasks)

        if n_tasks <= self.n_workers:
            # Few tasks, submit all at once
            return n_tasks

        # Aim for 2-4 chunks per worker
        optimal_chunks = self.n_workers * 3
        chunk_size = max(1, n_tasks // optimal_chunks)

        # Clamp to reasonable range
        chunk_size = max(5, min(chunk_size, 100))

        return chunk_size

    def _compute_scheduling_stats(
        self,
        original: List[ParallelTask],
        scheduled: List[ParallelTask]
    ) -> None:
        """Compute and store scheduling statistics."""
        # Check if order changed
        reordered = [t.task_id for t in original] != [t.task_id for t in scheduled]

        # Compute priority distribution
        priorities = [t.priority for t in scheduled]
        costs = [t.estimated_cost for t in scheduled]

        self._scheduling_stats = {
            "strategy": self.strategy.name,
            "n_tasks": len(scheduled),
            "reordered": reordered,
            "priority_range": (min(priorities), max(priorities)) if priorities else (0, 0),
            "cost_range": (min(costs), max(costs)) if costs else (0, 0),
            "avg_priority": sum(priorities) / len(priorities) if priorities else 0,
            "avg_cost": sum(costs) / len(costs) if costs else 0,
        }

        logger.debug(f"Scheduling stats: {self._scheduling_stats}")


def assign_task_priorities(
    tasks: List[ParallelTask],
    strategy: str = "nbody_first"
) -> None:
    """
    Assign priorities to tasks based on strategy.

    This is a convenience function that modifies tasks in-place.

    Parameters
    ----------
    tasks : List[ParallelTask]
        Tasks to assign priorities (modified in-place)
    strategy : str
        Priority strategy:
        - 'nbody_first': Higher n-body levels get higher priority
        - 'cost_first': Higher cost tasks get higher priority
        - 'balanced': Balance between cost and n-body level

    Examples
    --------
    >>> assign_task_priorities(tasks, strategy="nbody_first")
    >>> # Now tasks have priorities based on n-body level
    """
    if not tasks:
        return

    if strategy == "nbody_first":
        for task in tasks:
            task.priority = task.nbody * 10

    elif strategy == "cost_first":
        # Normalize costs to 0-100 range
        costs = [t.estimated_cost for t in tasks]
        max_cost = max(costs) if costs else 1.0
        min_cost = min(costs) if costs else 0.0
        cost_range = max_cost - min_cost if max_cost > min_cost else 1.0

        for task in tasks:
            normalized_cost = (task.estimated_cost - min_cost) / cost_range
            task.priority = int(normalized_cost * 100)

    elif strategy == "balanced":
        # Combine cost and n-body level
        costs = [t.estimated_cost for t in tasks]
        max_cost = max(costs) if costs else 1.0
        min_cost = min(costs) if costs else 0.0
        cost_range = max_cost - min_cost if max_cost > min_cost else 1.0

        max_nbody = max(t.nbody for t in tasks)

        for task in tasks:
            normalized_cost = (task.estimated_cost - min_cost) / cost_range
            normalized_nbody = task.nbody / max_nbody if max_nbody > 0 else 0

            # 60% weight on n-body, 40% on cost
            task.priority = int((normalized_nbody * 60 + normalized_cost * 40))

    else:
        logger.warning(f"Unknown priority strategy '{strategy}', priorities unchanged")

    logger.debug(f"Assigned priorities using strategy '{strategy}'")


def estimate_load_balance(
    tasks: List[ParallelTask],
    n_workers: int
) -> Dict[str, float]:
    """
    Estimate load balance for given task distribution.

    Parameters
    ----------
    tasks : List[ParallelTask]
        Tasks to analyze
    n_workers : int
        Number of workers

    Returns
    -------
    Dict[str, float]
        Load balance statistics:
        - 'ideal_time': Ideal time if perfectly balanced
        - 'estimated_time': Estimated actual time
        - 'efficiency': Load balance efficiency (0-1)
        - 'imbalance': Load imbalance factor (1 = perfect, >1 = imbalanced)
    """
    if not tasks or n_workers <= 0:
        return {
            "ideal_time": 0.0,
            "estimated_time": 0.0,
            "efficiency": 1.0,
            "imbalance": 1.0
        }

    # Simulate greedy assignment to workers
    worker_loads = [0.0] * n_workers

    for task in tasks:
        # Assign to least loaded worker
        min_worker = min(range(n_workers), key=lambda i: worker_loads[i])
        worker_loads[min_worker] += task.estimated_cost

    total_cost = sum(t.estimated_cost for t in tasks)
    ideal_time = total_cost / n_workers
    estimated_time = max(worker_loads)

    efficiency = ideal_time / estimated_time if estimated_time > 0 else 1.0
    imbalance = estimated_time / ideal_time if ideal_time > 0 else 1.0

    return {
        "ideal_time": ideal_time,
        "estimated_time": estimated_time,
        "efficiency": efficiency,
        "imbalance": imbalance
    }
