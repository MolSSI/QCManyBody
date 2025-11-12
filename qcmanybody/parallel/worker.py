"""
Worker functions for executing individual tasks.

This module contains the core task execution logic used by all executors.
Functions here are designed to be serializable for use with multiprocessing.
"""

import time
import traceback
import logging
from typing import Optional

from .task import ParallelTask, TaskResult, TaskStatus
from .base import ExecutorConfig

logger = logging.getLogger(__name__)


def execute_single_task(task: ParallelTask, config: ExecutorConfig) -> TaskResult:
    """Execute a single quantum chemistry task.

    This is the main worker function that all executors use. It:
    1. Validates the task
    2. Executes the QC calculation using QCEngine
    3. Handles errors and timeouts
    4. Returns a TaskResult

    Parameters
    ----------
    task : ParallelTask
        Task to execute
    config : ExecutorConfig
        Executor configuration

    Returns
    -------
    TaskResult
        Result of task execution (success or failure)

    Notes
    -----
    This function is designed to be pickle-able for use with
    multiprocessing and MPI. It has no side effects except logging.
    """
    import os

    worker_id = f"worker_{os.getpid()}"
    start_time = time.time()

    logger.debug(f"[{worker_id}] Starting task {task.task_id}")

    try:
        # Import QCEngine here to avoid issues with forking
        try:
            import qcengine as qcng
        except ImportError:
            return TaskResult(
                task_id=task.task_id,
                success=False,
                status=TaskStatus.FAILED,
                error_type="ImportError",
                error_message="qcengine not installed. Install with: pip install qcengine",
                execution_time=time.time() - start_time,
                worker_id=worker_id,
            )

        # Execute QC calculation
        # Extract program from atomic_input
        if hasattr(task.atomic_input, "model"):
            program = getattr(task.atomic_input.model, "program", "psi4")
        else:
            program = "psi4"  # Default fallback

        logger.debug(f"[{worker_id}] Running {program} for {task.task_id}")

        # Call QCEngine
        atomic_result = qcng.compute(task.atomic_input, program)

        execution_time = time.time() - start_time

        # Check if computation was successful
        if atomic_result.success:
            logger.debug(f"[{worker_id}] Task {task.task_id} succeeded in {execution_time:.2f}s")
            return TaskResult(
                task_id=task.task_id,
                success=True,
                status=TaskStatus.COMPLETED,
                atomic_result=atomic_result,
                execution_time=execution_time,
                worker_id=worker_id,
            )
        else:
            logger.warning(f"[{worker_id}] Task {task.task_id} failed: {atomic_result.error}")
            return TaskResult(
                task_id=task.task_id,
                success=False,
                status=TaskStatus.FAILED,
                error_type="QCEngineError",
                error_message=str(atomic_result.error),
                execution_time=execution_time,
                worker_id=worker_id,
            )

    except TimeoutError as e:
        logger.error(f"[{worker_id}] Task {task.task_id} timed out")
        return TaskResult(
            task_id=task.task_id,
            success=False,
            status=TaskStatus.TIMEOUT,
            error_type="TimeoutError",
            error_message=f"Task exceeded {config.timeout_per_task}s timeout",
            execution_time=time.time() - start_time,
            worker_id=worker_id,
        )

    except Exception as e:
        logger.error(f"[{worker_id}] Task {task.task_id} raised exception: {e}")
        return TaskResult(
            task_id=task.task_id,
            success=False,
            status=TaskStatus.FAILED,
            error_type=type(e).__name__,
            error_message=str(e),
            error_traceback=traceback.format_exc(),
            execution_time=time.time() - start_time,
            worker_id=worker_id,
        )


def execute_task(task: ParallelTask) -> TaskResult:
    """Execute a single task with default configuration.

    This is a convenience wrapper around execute_single_task that uses
    default configuration. It's designed for use with concurrent.futures
    and MPI executors where passing config is inconvenient.

    Parameters
    ----------
    task : ParallelTask
        Task to execute

    Returns
    -------
    TaskResult
        Result of task execution

    Notes
    -----
    This function uses default ExecutorConfig values. For more control,
    use execute_single_task directly with a custom config.
    """
    # Use default config
    config = ExecutorConfig()
    return execute_single_task(task, config)


def validate_environment() -> bool:
    """Validate that worker environment is properly configured.

    Returns
    -------
    bool
        True if environment is valid, False otherwise
    """
    try:
        import qcengine

        return True
    except ImportError:
        logger.error("qcengine not available in worker environment")
        return False
