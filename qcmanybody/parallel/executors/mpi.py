"""
MPI-based executor for distributed parallel execution.

This module provides an MPI-based executor for running many-body calculations
across multiple nodes in HPC clusters. It uses a master-worker architecture
where the master node distributes tasks and workers execute them.

NOTE: This is a skeleton implementation for Milestone 6. Full implementation
requires additional work on fault tolerance, load balancing, and optimization.

Requirements
------------
- mpi4py>=3.0
- MPI implementation (OpenMPI, MPICH, Intel MPI, etc.)

Installation
------------
pip install qcmanybody[mpi]

Usage
-----
mpirun -np 16 python your_script.py
"""

from typing import List, Optional, Callable
import logging
import os

from ..base import BaseParallelExecutor, ExecutorConfig
from ..task import ParallelTask, TaskResult

logger = logging.getLogger(__name__)

# Try to import mpi4py
try:
    from mpi4py import MPI
    MPI_AVAILABLE = True
except ImportError:
    MPI_AVAILABLE = False
    MPI = None


class MPIExecutor(BaseParallelExecutor):
    """
    MPI-based executor for distributed parallel execution.

    Uses a master-worker architecture:
    - **Master (rank 0)**: Distributes tasks and collects results
    - **Workers (rank 1-N)**: Execute tasks and return results

    **Master-Worker Communication Protocol:**

    1. Master sends task to worker
    2. Worker executes task
    3. Worker sends result back to master
    4. Repeat until all tasks complete

    **Advantages:**

    - Scales to hundreds/thousands of cores
    - Efficient inter-node communication
    - Works on HPC clusters with batch systems

    **Requirements:**

    - MPI implementation (OpenMPI, MPICH, etc.)
    - mpi4py Python package
    - Batch system configuration (SLURM, PBS, etc.)

    Parameters
    ----------
    config : ExecutorConfig
        Executor configuration

    Raises
    ------
    ImportError
        If mpi4py is not installed

    Examples
    --------
    >>> # In your Python script:
    >>> from qcmanybody.parallel.executors import MPIExecutor
    >>> from qcmanybody.parallel import ExecutorConfig
    >>>
    >>> config = ExecutorConfig(timeout_per_task=3600.0)
    >>> executor = MPIExecutor(config)
    >>>
    >>> with executor:
    ...     results = executor.execute(tasks)
    >>>
    >>> # Run with MPI:
    >>> # mpirun -np 16 python script.py

    Notes
    -----
    - This is a skeleton implementation (Milestone 6)
    - Full implementation requires additional work on:
      - Dynamic load balancing
      - Fault tolerance (node failures)
      - Checkpoint integration
      - Performance profiling
    """

    def __init__(self, config: Optional[ExecutorConfig] = None):
        """
        Initialize MPI executor.

        Parameters
        ----------
        config : ExecutorConfig, optional
            Executor configuration

        Raises
        ------
        ImportError
            If mpi4py is not installed
        RuntimeError
            If not running under MPI
        """
        if not MPI_AVAILABLE:
            raise ImportError(
                "mpi4py is required for MPIExecutor but is not installed.\n\n"
                "To use MPI-based parallel execution:\n"
                "1. Install a system MPI implementation (OpenMPI, MPICH, Intel MPI, etc.)\n"
                "   - Ubuntu/Debian: sudo apt-get install libopenmpi-dev openmpi-bin\n"
                "   - CentOS/RHEL: sudo yum install openmpi-devel\n"
                "   - macOS: brew install open-mpi\n"
                "2. Install qcmanybody with MPI support:\n"
                "   pip install qcmanybody[mpi]\n\n"
                "For more information, see the parallel execution documentation:\n"
                "https://molssi.github.io/QCManyBody/parallel_execution_guide.html"
            )

        super().__init__(config or ExecutorConfig())

        # Initialize MPI communicator
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        if self.size < 2:
            raise RuntimeError(
                "MPIExecutor requires at least 2 MPI processes (1 master + 1 worker).\n\n"
                "Current MPI size: 1 process\n\n"
                "To run with MPI:\n"
                "  mpirun -np <N> python your_script.py\n"
                "  # or with mpiexec:\n"
                "  mpiexec -n <N> python your_script.py\n\n"
                "Where <N> is the number of processes (e.g., 4, 8, 16).\n"
                "Example: mpirun -np 4 python your_script.py\n\n"
                "For HPC clusters, use your batch system (SLURM, PBS, etc.):\n"
                "  sbatch job_script.sh  # SLURM\n"
                "  qsub job_script.pbs   # PBS/Torque"
            )

        self.is_master = (self.rank == 0)
        self.n_workers = self.size - 1  # Exclude master

        if self.is_master:
            logger.info(
                f"MPI Master initialized: {self.size} processes "
                f"({self.n_workers} workers)"
            )
        else:
            logger.debug(f"MPI Worker {self.rank} initialized")

    def _initialize_resources(self) -> None:
        """
        Initialize MPI resources.

        For master: Set up task queue management
        For workers: Prepare to receive tasks
        """
        if self.is_master:
            logger.info("Master: Initializing resources")
            # Task queue management will be handled in _execute_tasks_impl
        else:
            logger.debug(f"Worker {self.rank}: Ready to receive tasks")

    def _cleanup_resources(self) -> None:
        """
        Clean up MPI resources.

        Ensures all workers are notified to shut down.
        """
        if self.is_master:
            logger.info("Master: Cleaning up resources")
            # Send shutdown signal to all workers
            for worker_rank in range(1, self.size):
                self.comm.send(None, dest=worker_rank, tag=0)
        else:
            logger.debug(f"Worker {self.rank}: Shutting down")

        # Barrier to ensure all processes reach cleanup
        self.comm.Barrier()

    def execute(
        self,
        tasks: List[ParallelTask],
        progress_callback: Optional[Callable[[str, int, int], None]] = None
    ) -> List[TaskResult]:
        """
        Execute tasks using MPI master-worker pattern.

        Master distributes tasks to workers and collects results.
        Workers execute tasks and return results.

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to execute
        progress_callback : Optional[Callable[[str, int, int], None]]
            Optional callback function called as:
            progress_callback(task_id, completed_count, total_count)

        Returns
        -------
        List[TaskResult]
            Results from task execution
        """
        if self.is_master:
            return self._master_execute(tasks, progress_callback)
        else:
            self._worker_loop()
            return []  # Workers don't return results

    def _master_execute(
        self,
        tasks: List[ParallelTask],
        progress_callback: Optional[Callable[[str, int, int], None]] = None
    ) -> List[TaskResult]:
        """
        Master process: Distribute tasks and collect results.

        Implements simple work queue pattern:
        1. Send initial tasks to all workers
        2. As workers complete, send more tasks
        3. Collect all results

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to execute

        Returns
        -------
        List[TaskResult]
            Collected results from all workers
        """
        if not tasks:
            return []

        logger.info(f"Master: Distributing {len(tasks)} tasks to {self.n_workers} workers")

        results = []
        task_queue = list(tasks)
        active_workers = set()

        # Send initial tasks to all workers
        for worker_rank in range(1, min(self.size, len(tasks) + 1)):
            if task_queue:
                task = task_queue.pop(0)
                logger.debug(f"Master: Sending task {task.task_id} to worker {worker_rank}")
                self.comm.send(task, dest=worker_rank, tag=0)
                active_workers.add(worker_rank)

        # Collect results and send more tasks
        while active_workers:
            # Receive result from any worker
            status = MPI.Status()
            result = self.comm.recv(source=MPI.ANY_SOURCE, tag=1, status=status)
            worker_rank = status.Get_source()

            results.append(result)
            logger.debug(
                f"Master: Received result for task {result.task_id} "
                f"from worker {worker_rank} "
                f"(success: {result.success})"
            )

            # Call progress callback
            if progress_callback:
                try:
                    progress_callback(result.task_id, len(results), len(tasks))
                except Exception as e:
                    logger.warning(f"Progress callback failed: {e}")

            # Send next task to this worker or mark idle
            if task_queue:
                task = task_queue.pop(0)
                logger.debug(f"Master: Sending task {task.task_id} to worker {worker_rank}")
                self.comm.send(task, dest=worker_rank, tag=0)
            else:
                # No more tasks, worker becomes idle
                active_workers.remove(worker_rank)
                logger.debug(f"Master: Worker {worker_rank} idle")

        logger.info(f"Master: Collected {len(results)} results")

        # Count successes/failures
        successes = sum(1 for r in results if r.success)
        failures = len(results) - successes
        logger.info(f"Master: {successes} succeeded, {failures} failed")

        return results

    def _worker_loop(self) -> None:
        """
        Worker process: Execute tasks until shutdown signal.

        Workers repeatedly:
        1. Receive task from master
        2. Execute task
        3. Send result back to master
        4. Repeat until receiving shutdown signal (None)
        """
        from ..worker import execute_task

        logger.debug(f"Worker {self.rank}: Entering work loop")

        while True:
            # Receive task from master
            task = self.comm.recv(source=0, tag=0)

            # None signals shutdown
            if task is None:
                logger.debug(f"Worker {self.rank}: Received shutdown signal")
                break

            logger.debug(f"Worker {self.rank}: Executing task {task.task_id}")

            # Execute task
            try:
                result = execute_task(task)
            except Exception as e:
                logger.error(f"Worker {self.rank}: Task {task.task_id} failed: {e}")
                result = TaskResult(
                    task_id=task.task_id,
                    success=False,
                    error=str(e),
                    result_data=None,
                    worker_id=f"mpi_rank_{self.rank}"
                )

            # Send result back to master
            logger.debug(f"Worker {self.rank}: Sending result for task {task.task_id}")
            self.comm.send(result, dest=0, tag=1)

        logger.debug(f"Worker {self.rank}: Exiting work loop")

    def get_worker_count(self) -> int:
        """
        Get number of MPI workers.

        Returns
        -------
        int
            Number of worker processes (size - 1)
        """
        return self.n_workers

    def is_initialized(self) -> bool:
        """
        Check if executor is initialized.

        Returns
        -------
        bool
            Always True for MPI (initialized in __init__)
        """
        return True


def create_mpi_executor(
    timeout_per_task: float = 3600.0,
    max_retries: int = 2
) -> MPIExecutor:
    """
    Convenience function to create an MPI executor.

    Parameters
    ----------
    timeout_per_task : float
        Timeout per task in seconds
    max_retries : int
        Maximum number of retries for failed tasks

    Returns
    -------
    MPIExecutor
        Configured MPI executor

    Examples
    --------
    >>> # Create MPI executor with custom timeout
    >>> executor = create_mpi_executor(timeout_per_task=7200.0)
    >>>
    >>> # Run with MPI:
    >>> # mpirun -np 16 python script.py
    """
    config = ExecutorConfig(
        timeout_per_task=timeout_per_task,
        max_retries=max_retries
    )
    return MPIExecutor(config)


# Utility function for checking MPI availability
def is_mpi_available() -> bool:
    """
    Check if MPI is available.

    Returns
    -------
    bool
        True if mpi4py is installed and we're running under MPI
    """
    if not MPI_AVAILABLE:
        return False

    try:
        comm = MPI.COMM_WORLD
        return comm.Get_size() > 1
    except Exception:
        return False


# Utility function for getting MPI info
def get_mpi_info() -> dict:
    """
    Get MPI configuration information.

    Returns
    -------
    dict
        MPI info including rank, size, processor name

    Raises
    ------
    RuntimeError
        If MPI is not available
    """
    if not MPI_AVAILABLE:
        raise RuntimeError("mpi4py not available")

    comm = MPI.COMM_WORLD
    return {
        "rank": comm.Get_rank(),
        "size": comm.Get_size(),
        "processor_name": MPI.Get_processor_name(),
        "is_master": comm.Get_rank() == 0,
        "n_workers": comm.Get_size() - 1
    }
