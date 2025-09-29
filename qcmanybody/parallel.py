"""Parallel execution engine for QCManyBody calculations.

This module implements level-by-level parallel execution of N-body calculations,
respecting mathematical dependencies while enabling parallelization within each level.

The parallel executor uses the dependency graph foundation from P1-002 to ensure
mathematical correctness while achieving performance improvements through parallelization.
"""

from __future__ import annotations

import logging
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from qcelemental.models import AtomicInput, AtomicResult, Molecule

from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import ManyBodyInput

try:
    from qcmanybody.computer import ManyBodyComputer
except Exception:  # pragma: no cover - circular import guard for documentation builds
    ManyBodyComputer = None  # type: ignore[assignment]

try:
    import qcengine as qcng
    HAS_QCENGINE = True
except ImportError:
    qcng = None
    HAS_QCENGINE = False

logger = logging.getLogger(__name__)

# Shared context used by multiprocessing workers. Each worker process receives a
# serialized copy of the configuration during initialization to avoid pickling
# the full executor state for every submitted task.
_WORKER_CONTEXT: Dict[str, Any] = {}


@dataclass
class ParallelConfig:
    """Configuration for parallel execution parameters.

    Attributes
    ----------
    max_workers : int
        Maximum number of parallel workers (default: 4)
    execution_mode : str
        Execution mode: "multiprocessing", "threading", or "serial" (default: "multiprocessing")
    memory_limit_mb : int
        Memory limit per worker in MB (default: 1000)
    timeout_seconds : int
        Timeout for individual fragment calculations in seconds (default: 3600)
    qcengine_config : dict
        Additional configuration for QCEngine execution
    """
    max_workers: int = 4
    execution_mode: str = "multiprocessing"
    memory_limit_mb: int = 1000
    timeout_seconds: int = 3600
    qcengine_config: Dict[str, Any] = None
    use_qcengine: bool = True
    qc_program: Optional[str] = "psi4"
    basis_set: Optional[str] = "sto-3g"
    default_driver: str = "energy"

    def __post_init__(self):
        if self.qcengine_config is None:
            self.qcengine_config = {}

        if self.use_qcengine and not HAS_QCENGINE:
            raise ImportError(
                "QCEngine is required for parallel execution but not available. "
                "Install with: pip install qcengine"
            )

        if self.execution_mode not in ["multiprocessing", "threading", "serial"]:
            raise ValueError(f"Invalid execution_mode: {self.execution_mode}")

        if self.max_workers < 1:
            raise ValueError(f"max_workers must be >= 1, got {self.max_workers}")


@dataclass
class FragmentTask:
    """Container describing a single fragment computation."""

    level: int
    model_chemistry: str
    label: str
    atomic_input: AtomicInput
    program: str


class ParallelManyBodyExecutor:
    """Parallel execution engine for many-body calculations.

    This class implements level-by-level parallel execution that respects N-body
    dependencies while enabling parallelization within each level. It builds on
    the P1-002 dependency graph foundation to ensure mathematical correctness.

    Parameters
    ----------
    core : ManyBodyCore
        The ManyBodyCore instance containing the dependency graph and calculation setup
    config : ParallelConfig
        Configuration for parallel execution parameters
    """

    # Class-level lock for thread-safe QCEngine initialization
    _qcengine_init_lock = threading.Lock()

    def __init__(
        self,
        core: ManyBodyCore,
        config: ParallelConfig,
        *,
        driver: Optional[str] = None,
        specifications: Optional[Dict[str, Dict[str, Any]]] = None,
    ):
        """Initialize the parallel executor.

        Parameters
        ----------
        core : ManyBodyCore
            ManyBodyCore instance with P1-002 dependency graph foundation
        config : ParallelConfig
            Parallel execution configuration
        driver : str
            Many-body driver requested (energy, gradient, hessian)
        specifications : Dict[str, Dict[str, Any]]
            Mapping of model chemistry labels to program + atomic specification dicts
        """
        self.core = core
        self.config = config
        self.driver = driver or config.default_driver
        self.specifications = dict(specifications or {})
        self._dependency_graph = core.dependency_graph

        # Validation: ensure P1-002 foundation is available
        if not hasattr(core, 'iterate_molecules_by_level'):
            raise RuntimeError(
                "ManyBodyCore missing iterate_molecules_by_level() method. "
                "P1-002 dependency graph foundation required."
            )

        # Track execution statistics
        self.execution_stats = {
            "total_fragments": 0,
            "levels_executed": 0,
            "parallel_time": 0.0,
            "sequential_time_estimate": 0.0,
            "speedup_factor": 0.0
        }

        logger.info(
            f"ParallelManyBodyExecutor initialized with {config.max_workers} workers "
            f"in {config.execution_mode} mode"
        )
        self._tasks_by_level = self._build_fragment_tasks()

    @classmethod
    def from_manybodyinput(cls, input_model: ManyBodyInput, config: ParallelConfig) -> "ParallelManyBodyExecutor":
        """Build an executor directly from a ``ManyBodyInput`` model."""

        if ManyBodyComputer is None:
            raise RuntimeError("ManyBodyComputer is unavailable; cannot construct executor from ManyBodyInput.")

        computer_model = ManyBodyComputer.from_manybodyinput(input_model, build_tasks=False)
        driver = computer_model.driver.value if hasattr(computer_model.driver, "value") else computer_model.driver
        specifications = cls._specifications_from_manybodyinput(input_model, driver)

        return cls(computer_model.qcmb_core, config, driver=driver, specifications=specifications)

    @staticmethod
    def _specifications_from_manybodyinput(
        input_model: ManyBodyInput, driver: str
    ) -> Dict[str, Dict[str, Any]]:
        """Convert ``ManyBodyInput`` atomic specifications into executor mapping."""

        specification_map: Dict[str, Dict[str, Any]] = {}

        for label, spec in input_model.specification.specification.items():
            spec_dict = spec.dict()
            program = spec_dict.pop("program")
            spec_dict["driver"] = driver
            spec_dict.pop("schema_name", None)
            spec_dict.pop("schema_version", None)

            specification_map[label] = {
                "program": program,
                "specification": spec_dict,
            }

        return specification_map

    def _build_fragment_tasks(self) -> Dict[int, List[FragmentTask]]:
        """Construct fragment tasks grouped by dependency level."""

        tasks: Dict[int, List[FragmentTask]] = {}

        for level, model_chemistry, label, molecule in self.core.iterate_molecules_by_level():
            spec_entry = self._resolve_specification(model_chemistry)
            spec_payload = deepcopy(spec_entry["specification"])

            # Ensure driver alignment
            spec_payload["driver"] = spec_payload.get("driver", self.driver)

            # Embed fragment-specific molecule
            atomic_input = AtomicInput(molecule=molecule, **spec_payload)

            # Handle embedding charges for supported programs
            if molecule.extras.get("embedding_charges") and spec_entry["program"] == "psi4":
                charges = molecule.extras["embedding_charges"]
                keywords = dict(atomic_input.keywords or {})
                fkw = dict(keywords.get("function_kwargs", {}))
                fkw.update({"external_potentials": charges})
                keywords["function_kwargs"] = fkw
                atomic_input = atomic_input.copy(update={"keywords": keywords})

            fragment_task = FragmentTask(
                level=level,
                model_chemistry=model_chemistry,
                label=label,
                atomic_input=atomic_input,
                program=spec_entry["program"],
            )

            tasks.setdefault(level, []).append(fragment_task)

        return tasks

    def _resolve_specification(self, model_chemistry: str) -> Dict[str, Any]:
        """Fetch the atomic specification entry for a model chemistry.

        Falls back to configuration defaults when an explicit specification
        is not provided (primarily for placeholder/demo workflows).
        """

        if model_chemistry in self.specifications:
            return self.specifications[model_chemistry]

        if not self.config.qc_program or not self.config.basis_set:
            raise KeyError(
                f"No specification available for model chemistry '{model_chemistry}' and no defaults configured."
            )

        # Build a minimal specification from config for backwards compatibility
        spec = {
            "program": self.config.qc_program,
            "specification": {
                "driver": self.driver or self.config.default_driver,
                "model": {"method": model_chemistry, "basis": self.config.basis_set},
                "keywords": {},
                "protocols": {},
                "extras": {},
            },
        }

        # Cache for future lookups to avoid recreating objects
        self.specifications[model_chemistry] = spec
        return spec

    def execute_fragment(self, task: FragmentTask) -> Tuple[str, AtomicResult]:
        """Execute a single fragment calculation using executor configuration."""

        return self._execute_fragment_static(task, self.config)

    @staticmethod
    def _execute_fragment_static(task: FragmentTask, config: ParallelConfig) -> Tuple[str, AtomicResult]:
        """Stateless fragment execution helper shared across worker contexts."""

        level = task.level
        model_chemistry = task.model_chemistry
        label = task.label
        atomic_input = task.atomic_input
        molecule = atomic_input.molecule

        try:
            logger.debug(f"Executing fragment {label} at level {level} with {model_chemistry}")

            if config.use_qcengine and HAS_QCENGINE:
                # Ensure QCEngine thread safety before execution
                ParallelManyBodyExecutor._ensure_qcengine_thread_safety()

                task_config_defaults = {
                    "memory": config.memory_limit_mb / 1024.0,
                    "ncores": 1,
                    "nnodes": 1,
                    "cores_per_rank": 1,
                    "retries": 0,
                }

                task_config = {
                    **task_config_defaults,
                    **config.qcengine_config.get("task_config", {}),
                }

                result = qcng.compute(
                    atomic_input,
                    task.program,
                    raise_error=True,
                    task_config=task_config
                )

                if not result.success:
                    raise RuntimeError(f"QCEngine calculation failed: {result.error}")

                logger.debug(f"Fragment {label} completed with energy {result.return_result}")

            else:
                logger.debug(f"Using placeholder calculation for {label} (QCEngine disabled)")

                time.sleep(0.01)

                natoms = len(molecule.symbols)
                energy_estimate = -natoms * 1.0 - sum(ord(c) for c in label) * 1e-6

                return_payload: Any
                properties: Dict[str, Any] = {"return_energy": energy_estimate}

                if atomic_input.driver == "energy":
                    return_payload = energy_estimate
                elif atomic_input.driver == "gradient":
                    gradient = np.zeros((natoms, 3))

                    return_payload = gradient
                    properties["return_gradient"] = gradient
                    properties["calcinfo_natom"] = natoms
                elif atomic_input.driver == "hessian":
                    gradient = np.zeros((natoms, 3))
                    size = 3 * natoms
                    hessian = np.zeros((size, size))
                    return_payload = hessian
                    properties["return_gradient"] = gradient
                    properties["return_hessian"] = hessian
                    properties["calcinfo_natom"] = natoms
                else:
                    raise RuntimeError(f"Unsupported driver for placeholder execution: {atomic_input.driver}")

                result = AtomicResult(
                    driver=atomic_input.driver,
                    model=atomic_input.model,
                    molecule=molecule,
                    return_result=return_payload,
                    success=True,
                    properties=properties,
                    provenance={"creator": "qcmanybody-parallel-placeholder", "version": "dev"}
                )

            logger.debug(f"Fragment {label} completed successfully")
            return label, result

        except Exception as exc:
            logger.error(f"Fragment {label} execution failed: {exc}")
            raise RuntimeError(f"Fragment calculation failed for {label}: {exc}")

    def execute_level_parallel(self, level: int, fragments_at_level: List[FragmentTask]) -> Dict[str, AtomicResult]:
        """Execute all fragments at a given level in parallel.

        This method executes all fragments at the same dependency level in parallel,
        ensuring mathematical correctness while maximizing performance.

        Parameters
        ----------
        level : int
            The N-body dependency level (1 for monomers, 2 for dimers, etc.)
        fragments_at_level : List[FragmentTask]
            List of fragment tasks at this level

        Returns
        -------
        Dict[str, AtomicResult]
            Dictionary mapping fragment labels to their calculation results

        Raises
        ------
        RuntimeError
            If any fragment calculation fails or parallel execution encounters errors
        """
        logger.info(f"Executing level {level} with {len(fragments_at_level)} fragments in parallel")

        if not fragments_at_level:
            logger.warning(f"No fragments found at level {level}")
            return {}

        level_results = {}
        start_time = time.time()

        if self.config.execution_mode == "serial":
            # Serial execution for debugging/comparison
            for fragment_task in fragments_at_level:
                label, result = self.execute_fragment(fragment_task)
                level_results[label] = result

        elif self.config.execution_mode == "threading":
            # Thread-based parallel execution
            with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                future_to_label = {
                    executor.submit(self.execute_fragment, fragment_task): fragment_task.label
                    for fragment_task in fragments_at_level
                }

                for future in as_completed(future_to_label, timeout=self.config.timeout_seconds):
                    label = future_to_label[future]
                    try:
                        result_label, result = future.result()
                        level_results[result_label] = result
                    except Exception as e:
                        logger.error(f"Fragment {label} failed in parallel execution: {e}")
                        raise RuntimeError(f"Parallel execution failed for fragment {label}: {e}")

        elif self.config.execution_mode == "multiprocessing":
            # Process-based parallel execution with spawn-safe worker initializer
            worker_config = deepcopy(self.config)

            with ProcessPoolExecutor(
                max_workers=self.config.max_workers,
                initializer=_multiprocessing_worker_initializer,
                initargs=(worker_config,)
            ) as executor:
                future_to_label = {
                    executor.submit(_execute_fragment_worker, fragment_task): fragment_task.label
                    for fragment_task in fragments_at_level
                }

                for future in as_completed(future_to_label, timeout=self.config.timeout_seconds):
                    label = future_to_label[future]
                    try:
                        result_label, result = future.result()
                        level_results[result_label] = result
                    except Exception as exc:
                        logger.error(f"Fragment {label} failed in parallel execution: {exc}")
                        raise RuntimeError(f"Parallel execution failed for fragment {label}: {exc}")

        execution_time = time.time() - start_time
        logger.info(f"Level {level} completed in {execution_time:.2f}s with {len(level_results)} results")

        # Update execution statistics
        self.execution_stats["total_fragments"] += len(fragments_at_level)
        self.execution_stats["parallel_time"] += execution_time

        return level_results

    def execute_full_calculation(self) -> Dict[str, AtomicResult]:
        """Execute complete many-body calculation with level-by-level parallelism.

        This method orchestrates the full parallel execution using the P1-002
        dependency graph foundation to ensure mathematical correctness.

        Returns
        -------
        Dict[str, AtomicResult]
            Dictionary mapping all fragment labels to their calculation results

        Raises
        ------
        RuntimeError
            If parallel execution fails or produces invalid results
        """
        logger.info("Starting full parallel many-body calculation")
        start_time = time.time()

        all_results: Dict[str, AtomicResult] = {}

        for level in sorted(self._tasks_by_level.keys()):
            fragments_at_level = self._tasks_by_level[level]
            level_results = self.execute_level_parallel(level, fragments_at_level)
            all_results.update(level_results)
            self.execution_stats["levels_executed"] += 1

        total_time = time.time() - start_time
        self.execution_stats["parallel_time"] = total_time

        # Estimate sequential time for speedup calculation
        # This is a rough estimate based on fragment count and average time per fragment
        avg_fragment_time = total_time / max(self.execution_stats["total_fragments"], 1)
        self.execution_stats["sequential_time_estimate"] = (
            avg_fragment_time * self.execution_stats["total_fragments"]
        )

        if self.execution_stats["parallel_time"] > 0:
            self.execution_stats["speedup_factor"] = (
                self.execution_stats["sequential_time_estimate"] / self.execution_stats["parallel_time"]
            )

        logger.info(f"Parallel calculation completed in {total_time:.2f}s")
        logger.info(f"Executed {self.execution_stats['total_fragments']} fragments "
                   f"across {self.execution_stats['levels_executed']} levels")
        logger.info(f"Estimated speedup: {self.execution_stats['speedup_factor']:.2f}x")

        return all_results

    def get_execution_statistics(self) -> Dict[str, Union[int, float]]:
        """Get detailed execution statistics.

        Returns
        -------
        Dict[str, Union[int, float]]
            Dictionary containing execution performance metrics
        """
        return self.execution_stats

    @staticmethod
    def _ensure_qcengine_thread_safety():
        """
        Ensure QCEngine global configuration is available in current thread.

        This method addresses the race condition in QCEngine's global configuration
        system when used in multithreaded environments. It uses a thread-safe
        approach to initialize QCEngine configuration if needed.
        """
        if not HAS_QCENGINE:
            return

        try:
            # Import qcengine in thread-safe manner
            import qcengine as qcng

            # Use class-level lock to prevent race conditions during initialization
            with ParallelManyBodyExecutor._qcengine_init_lock:
                try:
                    # Attempt to get configuration - this may trigger initialization
                    _ = qcng.get_config()
                    # Success - configuration is available
                    return
                except (KeyError, AttributeError) as e:
                    # Configuration missing - try to reinitialize
                    logger.debug(f"QCEngine config missing in thread, reinitializing: {e}")

                    try:
                        # Force QCEngine to reinitialize its global state
                        # This is a defensive approach for threading compatibility
                        from qcengine import config as qcng_config

                        # Try to trigger configuration reload
                        if hasattr(qcng_config, 'read_qcengine_task_environment'):
                            qcng_config.read_qcengine_task_environment()

                        # Verify configuration is now available
                        test_config = qcng.get_config()
                        logger.debug(f"QCEngine reinitialization successful: {test_config}")

                    except Exception as init_error:
                        # If reinitialization fails, we'll rely on task_config parameters
                        logger.warning(f"QCEngine thread initialization failed: {init_error}")
                        logger.warning("Will rely on explicit task_config parameters")

        except Exception as e:
            # If all else fails, log warning but continue - task_config should handle it
            logger.warning(f"QCEngine thread safety setup failed: {e}")
            logger.warning("Proceeding with task_config - may experience issues in threading mode")

    def validate_parallel_correctness(self, parallel_results: Dict[str, AtomicResult],
                                    sequential_results: Dict[str, AtomicResult],
                                    tolerance: float = 1e-12) -> bool:
        """Validate parallel results against sequential results.

        This method implements ultra-strict validation to ensure mathematical
        correctness of parallel execution results.

        Parameters
        ----------
        parallel_results : Dict[str, AtomicResult]
            Results from parallel execution
        sequential_results : Dict[str, AtomicResult]
            Results from sequential execution for comparison
        tolerance : float
            Numerical tolerance for comparison (default: 1e-12)

        Returns
        -------
        bool
            True if parallel results are within tolerance of sequential results

        Raises
        ------
        ValueError
            If results don't match within tolerance or have structural differences
        """
        logger.info(f"Validating parallel correctness with tolerance {tolerance}")

        # Check that we have the same number of results
        if len(parallel_results) != len(sequential_results):
            raise ValueError(
                f"Result count mismatch: parallel={len(parallel_results)}, "
                f"sequential={len(sequential_results)}"
            )

        # Check that all labels match
        parallel_labels = set(parallel_results.keys())
        sequential_labels = set(sequential_results.keys())
        if parallel_labels != sequential_labels:
            missing_parallel = sequential_labels - parallel_labels
            missing_sequential = parallel_labels - sequential_labels
            raise ValueError(
                f"Label mismatch: missing in parallel={missing_parallel}, "
                f"missing in sequential={missing_sequential}"
            )

        # Validate numerical results for each fragment
        max_difference = 0.0
        for label in parallel_results:
            parallel_result = parallel_results[label]
            sequential_result = sequential_results[label]

            # Compare return_result (energy values)
            if hasattr(parallel_result, 'return_result') and hasattr(sequential_result, 'return_result'):
                diff = abs(parallel_result.return_result - sequential_result.return_result)
                max_difference = max(max_difference, diff)

                if diff > tolerance:
                    raise ValueError(
                        f"Energy difference for {label} exceeds tolerance: "
                        f"{diff} > {tolerance}"
                    )

        logger.info(f"Parallel validation passed: max difference = {max_difference}")
        return True


def _multiprocessing_worker_initializer(config: ParallelConfig) -> None:
    """Initializer executed in each worker process to cache configuration."""

    _WORKER_CONTEXT["config"] = config


def _execute_fragment_worker(task: FragmentTask) -> Tuple[str, AtomicResult]:
    """ProcessPoolExecutor entry point that executes a fragment using cached config."""

    config = _WORKER_CONTEXT.get("config")
    if config is None:
        raise RuntimeError("Parallel worker received task before configuration initialization.")

    return ParallelManyBodyExecutor._execute_fragment_static(task, config)