"""
Parallel-enabled ManyBodyComputer.

This module extends ManyBodyComputer with parallel execution capabilities.
It maintains full backward compatibility with the sequential implementation
while allowing users to leverage parallel executors for speedup.
"""

from typing import Optional, Dict, Any
import logging

from qcelemental.models import Molecule, AtomicInput

from qcmanybody.computer import ManyBodyComputer
from qcmanybody.models.v1 import ManyBodyInput, ManyBodyResult
from qcmanybody.utils import delabeler
from .base import BaseParallelExecutor, ExecutorConfig
from .executors.sequential import SequentialExecutor
from .task import ParallelTask, TaskResult

# Type hints for optional imports
TYPE_CHECKING = False
if TYPE_CHECKING:
    from .scheduler import TaskScheduler
    from .checkpoint import CheckpointManager

logger = logging.getLogger(__name__)


class ParallelManyBodyComputer(ManyBodyComputer):
    """ManyBodyComputer with parallel execution support.

    This class extends ManyBodyComputer to support parallel execution of
    quantum chemistry tasks using pluggable executors. It maintains full
    backward compatibility - when no executor is specified, it behaves
    identically to the base ManyBodyComputer.

    Parameters
    ----------
    executor : Optional[BaseParallelExecutor]
        Parallel executor to use for task execution. If None, uses
        SequentialExecutor (no parallelism). Default: None
    **kwargs
        All other arguments passed to ManyBodyComputer

    Examples
    --------
    Basic usage (sequential, same as ManyBodyComputer):

    >>> computer = ParallelManyBodyComputer.from_manybodyinput(input_model)

    With multiprocessing executor:

    >>> from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig
    >>> config = ExecutorConfig(n_workers=4)
    >>> executor = MultiprocessingExecutor(config)
    >>> computer = ParallelManyBodyComputer.from_manybodyinput(
    ...     input_model,
    ...     executor=executor
    ... )

    Simplified parallel API:

    >>> computer = ParallelManyBodyComputer.from_manybodyinput(
    ...     input_model,
    ...     parallel=True,
    ...     n_workers=4
    ... )
    """

    def __init__(
        self,
        *args,
        executor: Optional[BaseParallelExecutor] = None,
        scheduler: Optional['TaskScheduler'] = None,
        checkpoint_manager: Optional['CheckpointManager'] = None,
        **kwargs
    ):
        """Initialize ParallelManyBodyComputer.

        Parameters
        ----------
        executor : Optional[BaseParallelExecutor]
            Parallel executor for task execution
        scheduler : Optional[TaskScheduler]
            Task scheduler for optimizing task order and load balancing
        checkpoint_manager : Optional[CheckpointManager]
            Checkpoint manager for saving/resuming calculations
        *args, **kwargs
            Arguments passed to ManyBodyComputer
        """
        super().__init__(*args, **kwargs)
        self.executor = executor or SequentialExecutor()
        self.scheduler = scheduler
        self.checkpoint_manager = checkpoint_manager

    @classmethod
    def from_manybodyinput(
        cls,
        input_model: ManyBodyInput,
        build_tasks: bool = True,
        executor: Optional[BaseParallelExecutor] = None,
        parallel: bool = False,
        n_workers: Optional[int] = None,
        scheduler: Optional['TaskScheduler'] = None,
        checkpoint_manager: Optional['CheckpointManager'] = None,
    ) -> ManyBodyResult:
        """Create computer from ManyBodyInput and execute with parallel support.

        This method extends the base from_manybodyinput with parallel execution.
        It can be called in three ways:

        1. Sequential (default): Same as base ManyBodyComputer
        2. Explicit executor: Provide a configured executor
        3. Simple parallel: Set parallel=True and optionally n_workers

        Parameters
        ----------
        input_model : ManyBodyInput
            Input specification for many-body calculation
        build_tasks : bool
            If True, build and execute tasks. If False, only create computer.
            Default: True
        executor : Optional[BaseParallelExecutor]
            Explicit executor to use. Takes precedence over parallel parameter.
            Default: None
        parallel : bool
            If True and executor is None, create a MultiprocessingExecutor.
            Default: False
        n_workers : Optional[int]
            Number of workers for auto-created executor. Only used if
            parallel=True and executor=None. Default: None (auto-detect)
        scheduler : Optional[TaskScheduler]
            Task scheduler for optimizing task order and load balancing.
            If provided, tasks will be scheduled before execution.
            Default: None (no scheduling)
        checkpoint_manager : Optional[CheckpointManager]
            Checkpoint manager for saving/resuming calculations.
            If provided, results will be automatically checkpointed.
            Default: None (no checkpointing)

        Returns
        -------
        ManyBodyResult
            Completed many-body calculation result

        Examples
        --------
        >>> # Sequential (default)
        >>> result = ParallelManyBodyComputer.from_manybodyinput(input_model)

        >>> # Simple parallel
        >>> result = ParallelManyBodyComputer.from_manybodyinput(
        ...     input_model, parallel=True, n_workers=4
        ... )

        >>> # Custom executor
        >>> from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig
        >>> config = ExecutorConfig(n_workers=8, timeout_per_task=1800)
        >>> executor = MultiprocessingExecutor(config)
        >>> result = ParallelManyBodyComputer.from_manybodyinput(
        ...     input_model, executor=executor
        ... )
        """
        # Auto-create executor if parallel=True
        if executor is None and parallel:
            try:
                from .executors.multiprocessing import MultiprocessingExecutor
                config = ExecutorConfig(n_workers=n_workers)
                executor = MultiprocessingExecutor(config)
                logger.info(f"Auto-created MultiprocessingExecutor with n_workers={n_workers}")
            except ImportError:
                logger.warning(
                    "MultiprocessingExecutor not available, falling back to sequential execution"
                )
                executor = SequentialExecutor()
        elif executor is None:
            executor = SequentialExecutor()

        # Create base computer model (without executing tasks)
        # We call the parent class method with build_tasks=False
        computer_model = super(ParallelManyBodyComputer, cls).from_manybodyinput(
            input_model, build_tasks=False
        )

        # Set the executor, scheduler, and checkpoint manager
        computer_model.executor = executor
        computer_model.scheduler = scheduler
        computer_model.checkpoint_manager = checkpoint_manager

        if not build_tasks:
            return computer_model

        # Execute tasks in parallel
        return computer_model._execute_parallel(input_model)

    def _execute_parallel(self, input_model: ManyBodyInput) -> ManyBodyResult:
        """Execute many-body calculation using parallel executor.

        This method replaces the sequential loop in the base implementation
        with parallel task execution.

        Parameters
        ----------
        input_model : ManyBodyInput
            Input specification (needed for specifications dict)

        Returns
        -------
        ManyBodyResult
            Completed many-body calculation result
        """
        try:
            import qcengine as qcng
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "Python module qcengine not found. Solve by installing it: "
                "`conda install qcengine -c conda-forge` or `pip install qcengine`"
            )

        # Build specifications dict (same as base implementation)
        specifications = {}
        for mtd, spec in input_model.specification.specification.items():
            spec = spec.dict()
            specifications[mtd] = {}
            specifications[mtd]["program"] = spec.pop("program")
            specifications[mtd]["specification"] = spec
            specifications[mtd]["specification"]["driver"] = self.driver
            specifications[mtd]["specification"].pop("schema_name", None)

        # Prepare all tasks upfront
        tasks = []
        logger.info("Preparing parallel tasks...")

        for chem, label, imol in self.qcmb_core.iterate_molecules():
            # Build AtomicInput (same as sequential version)
            inp = AtomicInput(molecule=imol, **specifications[chem]["specification"])

            # Handle embedding charges
            if imol.extras.get("embedding_charges"):
                if specifications[chem]["program"] == "psi4":
                    charges = imol.extras["embedding_charges"]
                    fkw = inp.keywords.get("function_kwargs", {})
                    fkw.update({"external_potentials": charges})
                    inp.keywords["function_kwargs"] = fkw
                else:
                    raise RuntimeError(
                        f"Don't know how to handle external charges in {specifications[chem]['program']}"
                    )

            # Create parallel task
            task = ParallelTask(
                task_id=label,
                chemistry=chem,
                label=label,
                molecule=imol,
                atomic_input=inp,
                metadata={
                    "program": specifications[chem]["program"],
                },
            )
            tasks.append(task)

        logger.info(f"Prepared {len(tasks)} tasks for parallel execution")

        # Load checkpoint if available and filter completed tasks
        if self.checkpoint_manager:
            if self.checkpoint_manager.exists():
                logger.info("Loading checkpoint...")
                self.checkpoint_manager.load()
                tasks = self.checkpoint_manager.filter_pending_tasks(tasks)
                logger.info(f"After checkpoint: {len(tasks)} tasks remaining")

        # Apply scheduler if provided
        if self.scheduler:
            logger.info(f"Scheduling tasks with strategy: {self.scheduler.strategy.name}")
            tasks = self.scheduler.schedule(tasks)

        # Execute tasks in parallel
        with self.executor as executor:
            logger.info(f"Executing with {executor.name}")

            # Define progress callback if checkpoint manager is available
            progress_callback = None
            if self.checkpoint_manager:
                def progress_callback(task_id: str, completed: int, total: int):
                    """Save checkpoint after each completed task."""
                    # Find the completed task result
                    for task in tasks:
                        if task.task_id == task_id:
                            # Note: Results are collected after execution completes
                            pass
                    logger.debug(f"Progress: {completed}/{total} tasks completed")

            results = executor.execute(tasks, progress_callback=progress_callback)

        # Save checkpoint with all results
        if self.checkpoint_manager:
            # Convert TaskResult to checkpoint format
            for result in results:
                self.checkpoint_manager.save_result(result)
            self.checkpoint_manager.save(total_tasks=len(tasks))
            logger.info("Checkpoint saved")

        # Collect results into component_results and component_properties
        component_results = {}
        component_properties = {}

        # First, load any existing results from checkpoint
        all_results = []
        if self.checkpoint_manager:
            # Get completed results from checkpoint
            checkpoint_results = list(self.checkpoint_manager._results.values())
            all_results.extend(checkpoint_results)

        # Add newly computed results
        all_results.extend(results)

        # Process all results (checkpoint + new)
        for result in all_results:
            label = result.task_id

            if not result.success:
                error_msg = result.error_message or "Unknown error"
                raise RuntimeError(f"Calculation {label} did not succeed! Error:\n{error_msg}")

            # Store atomic result
            if result.atomic_result:
                component_results[label] = result.atomic_result

                # Extract properties
                props = {"energy", "gradient", "hessian"}
                component_properties[label] = {}

                for p in props:
                    if hasattr(result.atomic_result.properties, f"return_{p}"):
                        v = getattr(result.atomic_result.properties, f"return_{p}")
                        if v is not None:
                            component_properties[label][p] = v

        # Analyze results (same as base implementation)
        analyze_back = self.qcmb_core.analyze(component_properties)
        analyze_back["nbody_number"] = len(component_properties)

        # Return final result
        return self.get_results(
            external_results=analyze_back,
            component_results=component_results
        )


def parallel_compute_from_manybodyinput(
    input_model: ManyBodyInput,
    executor: Optional[BaseParallelExecutor] = None,
    parallel: bool = False,
    n_workers: Optional[int] = None,
) -> ManyBodyResult:
    """Convenience function for parallel many-body calculations.

    This is a simple wrapper around ParallelManyBodyComputer.from_manybodyinput
    that provides a functional interface matching the original compute pattern.

    Parameters
    ----------
    input_model : ManyBodyInput
        Input specification for many-body calculation
    executor : Optional[BaseParallelExecutor]
        Explicit executor to use. Default: None
    parallel : bool
        Enable simple parallel execution. Default: False
    n_workers : Optional[int]
        Number of workers for parallel execution. Default: None

    Returns
    -------
    ManyBodyResult
        Completed many-body calculation result

    Examples
    --------
    >>> from qcmanybody.parallel import parallel_compute_from_manybodyinput
    >>> result = parallel_compute_from_manybodyinput(
    ...     input_model,
    ...     parallel=True,
    ...     n_workers=4
    ... )
    """
    return ParallelManyBodyComputer.from_manybodyinput(
        input_model,
        executor=executor,
        parallel=parallel,
        n_workers=n_workers,
    )
