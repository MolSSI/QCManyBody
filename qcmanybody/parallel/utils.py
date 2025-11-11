"""
Utility functions for parallel execution.
"""

import logging
from typing import List, Dict, Any
from .task import ParallelTask, TaskResult

logger = logging.getLogger(__name__)


def estimate_task_cost(task: ParallelTask) -> float:
    """Estimate the computational cost of a task.

    This uses simple heuristics based on:
    - Number of atoms
    - QC method
    - Basis set size

    Parameters
    ----------
    task : ParallelTask
        Task to estimate cost for

    Returns
    -------
    float
        Estimated relative cost (higher = more expensive)
    """
    # Get number of atoms
    natoms = len(task.molecule.symbols)

    # Method scaling factors (very rough estimates)
    method_costs = {
        "hf": 1.0,
        "mp2": 5.0,
        "ccsd": 20.0,
        "ccsd(t)": 100.0,
        "b3lyp": 2.0,
        "wb97x": 3.0,
    }

    # Extract method from chemistry string
    method = task.chemistry.split("/")[0].lower()
    method_cost = method_costs.get(method, 10.0)  # Default moderate cost

    # Atom scaling (roughly N^3 to N^4 for most methods)
    atom_cost = natoms**3.5

    return method_cost * atom_cost


def assign_priorities(tasks: List[ParallelTask], strategy: str = "nbody") -> List[ParallelTask]:
    """Assign priorities to tasks based on a strategy.

    Parameters
    ----------
    tasks : List[ParallelTask]
        Tasks to assign priorities to
    strategy : str
        Priority strategy. Options:
        - "nbody": Lower n-body first (default)
        - "cost": Cheaper tasks first
        - "cost_desc": Expensive tasks first

    Returns
    -------
    List[ParallelTask]
        Tasks with updated priorities
    """
    if strategy == "nbody":
        # Lower n-body = higher priority (negative so higher nbody is lower priority)
        for task in tasks:
            task.priority = -task.nbody

    elif strategy == "cost":
        # Cheaper tasks first
        costs = [estimate_task_cost(t) for t in tasks]
        max_cost = max(costs) if costs else 1.0
        for task, cost in zip(tasks, costs):
            task.priority = int(100 * (1.0 - cost / max_cost))

    elif strategy == "cost_desc":
        # Expensive tasks first
        costs = [estimate_task_cost(t) for t in tasks]
        max_cost = max(costs) if costs else 1.0
        for task, cost in zip(tasks, costs):
            task.priority = int(100 * (cost / max_cost))

    else:
        raise ValueError(f"Unknown priority strategy: {strategy}")

    return tasks


def compute_execution_stats(results: List[TaskResult]) -> Dict[str, Any]:
    """Compute execution statistics from results.

    Parameters
    ----------
    results : List[TaskResult]
        Task results

    Returns
    -------
    Dict[str, Any]
        Statistics dictionary with:
        - total_tasks: Total number of tasks
        - successful: Number of successful tasks
        - failed: Number of failed tasks
        - total_time: Sum of all execution times
        - avg_time: Average execution time
        - min_time: Minimum execution time
        - max_time: Maximum execution time
        - success_rate: Fraction of successful tasks
    """
    if not results:
        return {
            "total_tasks": 0,
            "successful": 0,
            "failed": 0,
            "total_time": 0.0,
            "avg_time": 0.0,
            "min_time": 0.0,
            "max_time": 0.0,
            "success_rate": 0.0,
        }

    total_tasks = len(results)
    successful = sum(1 for r in results if r.success)
    failed = total_tasks - successful

    execution_times = [r.execution_time for r in results]
    total_time = sum(execution_times)
    avg_time = total_time / total_tasks if total_tasks > 0 else 0.0
    min_time = min(execution_times) if execution_times else 0.0
    max_time = max(execution_times) if execution_times else 0.0

    success_rate = successful / total_tasks if total_tasks > 0 else 0.0

    return {
        "total_tasks": total_tasks,
        "successful": successful,
        "failed": failed,
        "total_time": total_time,
        "avg_time": avg_time,
        "min_time": min_time,
        "max_time": max_time,
        "success_rate": success_rate,
    }


def format_execution_summary(results: List[TaskResult]) -> str:
    """Format execution statistics as a human-readable string.

    Parameters
    ----------
    results : List[TaskResult]
        Task results

    Returns
    -------
    str
        Formatted summary string
    """
    stats = compute_execution_stats(results)

    summary = f"""
Parallel Execution Summary:
---------------------------
Total tasks:    {stats['total_tasks']}
Successful:     {stats['successful']} ({stats['success_rate']*100:.1f}%)
Failed:         {stats['failed']}

Timing:
  Total time:   {stats['total_time']:.2f}s
  Average:      {stats['avg_time']:.2f}s
  Min:          {stats['min_time']:.2f}s
  Max:          {stats['max_time']:.2f}s
"""

    return summary.strip()
