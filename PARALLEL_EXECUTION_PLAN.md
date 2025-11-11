# QCManyBody Parallel Execution Module - Development Plan

**Version:** 1.0
**Date:** 2025-11-11
**Status:** Planning Phase
**Owner:** Development Team

---

## Executive Summary

This document outlines the comprehensive development plan for adding parallel execution capabilities to QCManyBody. The module will enable efficient computation of many-body expansion calculations by parallelizing independent quantum chemistry tasks across multiple cores and eventually across distributed systems using MPI.

**Current Bottleneck:** Sequential execution at `computer.py:465` - each QC calculation blocks until completion.

**Target:** 10-100x speedup for typical MBE calculations depending on system size and available resources.

**CLI Integration Status (2025-11-11):** The CLI has been merged into main. Parallel execution must now be integrated with the CLI to allow users to specify parallel execution options via input files and command-line arguments.

---

## Table of Contents

1. [Architecture Vision](#architecture-vision)
2. [Project Milestones](#project-milestones)
3. [Development Sprints](#development-sprints)
4. [Epic Breakdown](#epic-breakdown)
5. [Technical Specifications](#technical-specifications)
6. [Teaching & Documentation](#teaching--documentation)
7. [Testing Strategy](#testing-strategy)
8. [Risk Management](#risk-management)
9. [Success Metrics](#success-metrics)

---

## Architecture Vision

### Design Principles

1. **Non-invasive:** Existing sequential code remains fully functional
2. **Pluggable:** Parallel executors implement common interface
3. **Progressive:** Start with multiprocessing, extend to MPI
4. **Backward Compatible:** All existing APIs continue to work
5. **Type-safe:** Full type hints throughout

### Module Structure

```
qcmanybody/
├── parallel/
│   ├── __init__.py              # Public API exports
│   ├── base.py                  # BaseParallelExecutor abstract class
│   ├── executors/
│   │   ├── __init__.py
│   │   ├── sequential.py        # Default sequential (no-op wrapper)
│   │   ├── multiprocessing.py   # Python multiprocessing executor
│   │   ├── futures.py           # concurrent.futures executor
│   │   └── mpi.py               # MPI executor (future)
│   ├── task.py                  # ParallelTask, TaskResult classes
│   ├── scheduler.py             # Task scheduling & load balancing
│   ├── compute_parallel.py      # ParallelManyBodyComputer class
│   ├── core_parallel.py         # ParallelManyBodyCore extensions
│   └── utils.py                 # Parallel utilities
├── core.py                      # Extended with parallel support
├── computer.py                  # Extended with parallel support
└── builder.py                   # Enhanced for parallel task metadata
```

### Key Classes

```python
# Base executor interface
class BaseParallelExecutor(ABC):
    @abstractmethod
    def execute(self, tasks: List[ParallelTask]) -> List[TaskResult]: ...

    @abstractmethod
    def shutdown(self) -> None: ...

# Parallel task specification
@dataclass
class ParallelTask:
    task_id: str
    chemistry: str
    label: str
    molecule: Molecule
    specification: AtomicInput
    priority: int = 0

# Extended computer
class ParallelManyBodyComputer(ManyBodyComputer):
    def __init__(self, ..., executor: Optional[BaseParallelExecutor] = None):
        self.executor = executor or SequentialExecutor()
```

---

## Project Milestones

### Milestone 1: Foundation & Planning ✅ COMPLETE
**Duration:** 2 weeks
**Status:** Complete (2025-11-11)

**Deliverables:**
- ✅ Codebase analysis complete
- ✅ Development plan created
- ✅ Architecture design finalized
- ✅ Parallel module structure created
- ✅ CLI merged into main branch

**Exit Criteria:**
- ✅ All team members understand current architecture
- ✅ Parallel module structure agreed upon
- ✅ CLI integration needs identified

---

### Milestone 2: Skeleton & Infrastructure ✅ COMPLETE
**Duration:** 3 weeks
**Status:** Complete (2025-11-11)

**Deliverables:**
- ✅ `/qcmanybody/parallel/` folder structure created
- ✅ `BaseParallelExecutor` abstract class implemented
- ✅ `SequentialExecutor` (no-op wrapper) implemented
- ✅ `MultiprocessingExecutor` implemented
- ✅ `ParallelTask` and `TaskResult` data models
- ✅ `ParallelManyBodyComputer` class implemented
- ✅ Worker functions implemented
- [ ] Basic test harness for executors
- [ ] CI integration for parallel module

**Exit Criteria:**
- ✅ All base classes pass type checking
- ✅ Core parallel execution infrastructure complete
- [ ] 80%+ code coverage for base classes

---

### Milestone 3: CLI Integration & Testing 🎯 CURRENT
**Duration:** 2 weeks
**Target Date:** Week 7
**Status:** In Progress

**Deliverables:**
- [ ] Add `ExecutionSchema` to CLI input schema for parallel configuration
- [ ] Update `run.py` to use `ParallelManyBodyComputer`
- [ ] Add `--parallel` and `--n-workers` CLI arguments
- [ ] Add parallel execution examples (JSON/YAML input files)
- [ ] Basic test harness for executors
- [ ] Integration tests with CLI
- [ ] Update CLI documentation for parallel execution
- [ ] Performance benchmarks vs. sequential

**Exit Criteria:**
- CLI can execute calculations with parallel flag
- All existing CLI tests pass with parallel execution
- 2x+ speedup on 4-core system for CLI-based runs
- Documentation includes parallel execution examples

---

### Milestone 4: Multiprocessing Executor (Completion)
**Duration:** 2 weeks
**Target Date:** Week 9

**Deliverables:**
- ✅ `MultiprocessingExecutor` fully implemented
- ✅ Task serialization/deserialization
- [ ] Advanced error handling and task retry logic
- [ ] Progress tracking and logging enhancements
- [ ] Comprehensive unit tests
- [ ] Memory leak testing
- [ ] Error scenario testing

**Exit Criteria:**
- All existing tests pass with MultiprocessingExecutor
- No memory leaks in long-running calculations
- 90%+ code coverage for parallel module

---

### Milestone 5: Advanced Features & Optimization
**Duration:** 3 weeks
**Target Date:** Week 12

**Deliverables:**
- [ ] Dynamic load balancing
- [ ] Task priority scheduling
- [ ] Chunked task submission (memory optimization)
- [ ] Result caching and checkpointing
- [ ] `concurrent.futures` executor variant
- [ ] Performance profiling tools
- [ ] Optimization guide documentation

**Exit Criteria:**
- Load balancing reduces total time by 15%+
- Checkpoint/restart works reliably
- Memory usage scales linearly with worker count

---

### Milestone 6: MPI Support (Distributed Computing)
**Duration:** 4 weeks
**Target Date:** Week 16

**Deliverables:**
- [ ] `mpi4py` integration and dependency management
- [ ] `MPIExecutor` implementation
- [ ] Master-worker communication protocol
- [ ] Fault tolerance for node failures
- [ ] HPC cluster deployment guide
- [ ] SLURM/PBS job submission templates
- [ ] Multi-node performance benchmarks

**Exit Criteria:**
- Scales to 100+ workers across multiple nodes
- Graceful degradation on worker failure
- Documentation for 3+ HPC systems

---

### Milestone 7: Production Readiness
**Duration:** 2 weeks
**Target Date:** Week 18

**Deliverables:**
- [ ] Full API documentation
- [ ] User tutorials (beginner to advanced)
- [ ] Migration guide from sequential to parallel
- [ ] Performance tuning guide
- [ ] Security audit (pickle safety, etc.)
- [ ] Release notes and changelog
- [ ] PyPI package update

**Exit Criteria:**
- 95%+ code coverage across parallel module
- Zero P0/P1 bugs
- Documentation reviewed by 2+ external users
- Performance within 10% of theoretical maximum

---

## Development Sprints

### Sprint 1: Setup & Base Classes (Weeks 1-2)
**Goal:** Establish parallel module foundation

**Stories:**

1. **[SETUP-1] Create parallel module structure**
   - Create `/qcmanybody/parallel/` directory
   - Add `__init__.py` with module docstring
   - Create subdirectories: `executors/`, `tests/`
   - Update root `__init__.py` exports
   - **Points:** 2

2. **[BASE-1] Implement ParallelTask dataclass**
   - Define `ParallelTask` with all required fields
   - Add validation methods
   - Implement `__hash__` and `__eq__` for task deduplication
   - Add type hints and docstrings
   - **Points:** 3

3. **[BASE-2] Implement TaskResult dataclass**
   - Define `TaskResult` for success/failure states
   - Include timing and resource usage metadata
   - Add serialization helpers
   - **Points:** 2

4. **[BASE-3] Create BaseParallelExecutor ABC**
   - Define abstract interface
   - Add context manager support (`__enter__`/`__exit__`)
   - Document executor contract
   - **Points:** 3

5. **[TEST-1] Set up test infrastructure**
   - Create `conftest.py` with fixtures
   - Add mock QC program for testing
   - Create test molecule library
   - **Points:** 3

6. **[DOC-1] Write architecture documentation**
   - Document design decisions
   - Create class diagrams
   - Write contribution guide for parallel module
   - **Points:** 2

**Total Points:** 15

---

### Sprint 2: Sequential Executor & Integration (Weeks 3-4)
**Goal:** Implement no-op executor to validate interface

**Stories:**

1. **[SEQ-1] Implement SequentialExecutor**
   - Inherit from `BaseParallelExecutor`
   - Execute tasks one-by-one in current thread
   - Match exact behavior of original code
   - **Points:** 3

2. **[SEQ-2] Add comprehensive executor tests**
   - Test successful task execution
   - Test error handling
   - Test context manager behavior
   - Verify resource cleanup
   - **Points:** 5

3. **[INT-1] Extend ManyBodyCore for parallel**
   - Add `prepare_parallel_tasks()` method
   - Extract task generation logic from `iterate_molecules()`
   - Keep backward compatibility
   - **Points:** 5

4. **[INT-2] Create ParallelManyBodyComputer stub**
   - Inherit from `ManyBodyComputer`
   - Add `executor` parameter
   - Implement constructor
   - **Points:** 3

5. **[TEST-2] Validation against original**
   - Compare results bit-for-bit
   - Test all BSSE types (cp, nocp, vmfc)
   - Test all drivers (energy, gradient, hessian)
   - **Points:** 5

**Total Points:** 21

---

### Sprint 2.5: CLI Integration (Week 5) 🎯 CURRENT SPRINT
**Goal:** Integrate parallel execution with CLI

**Stories:**

1. **[CLI-1] Add ExecutionSchema to input schema**
   - Create `ExecutionSchema` class in `input_schema.py`
   - Add fields: `parallel`, `n_workers`, `executor_type`, `timeout`
   - Add to `QCManyBodyInput` model
   - Update schema validation
   - **Points:** 3

2. **[CLI-2] Update CLI run command**
   - Modify `run.py` to import ParallelManyBodyComputer
   - Extract execution config from CLI input
   - Pass execution parameters to computer
   - Handle both parallel and sequential modes
   - **Points:** 5

3. **[CLI-3] Add command-line arguments**
   - Add `--parallel` flag to run command
   - Add `--n-workers` argument
   - Add `--executor-type` argument
   - CLI args override input file settings
   - **Points:** 2

4. **[CLI-4] Create parallel execution examples**
   - Add `examples/cli/parallel_basic.json`
   - Add `examples/cli/parallel_custom.yaml`
   - Update examples README
   - **Points:** 2

5. **[CLI-5] Update CLI documentation**
   - Add parallel execution section to `cli_guide.md`
   - Document execution schema fields
   - Document CLI arguments
   - Add troubleshooting section
   - **Points:** 3

6. **[TEST-CLI] Integration tests**
   - Test CLI with parallel execution enabled
   - Test CLI argument overrides
   - Verify results match sequential
   - **Points:** 5

**Total Points:** 20

---

### Sprint 3: Multiprocessing Core (Weeks 5-6)
**Goal:** Implement basic multiprocessing executor

**Stories:**

1. **[MP-1] Design task serialization**
   - Handle molecule pickling
   - Handle QCEngine input serialization
   - Test with all molecule types
   - **Points:** 5

2. **[MP-2] Implement MultiprocessingExecutor**
   - Use `multiprocessing.Pool`
   - Implement worker function
   - Add process initialization
   - **Points:** 8

3. **[MP-3] Add basic error handling**
   - Catch and wrap worker exceptions
   - Timeout handling
   - Process crash detection
   - **Points:** 5

4. **[TEST-3] Multiprocessing correctness tests**
   - Verify results match sequential
   - Test with 2, 4, 8 workers
   - Test various molecular systems
   - **Points:** 5

**Total Points:** 23

---

### Sprint 4: Progress Tracking & Logging (Weeks 7-8)
**Goal:** Add visibility into parallel execution

**Stories:**

1. **[LOG-1] Design logging strategy**
   - Worker-safe logging configuration
   - Structured log format
   - Performance counters
   - **Points:** 3

2. **[LOG-2] Implement progress tracking**
   - Add `ProgressTracker` class
   - Integrate with executors
   - Support callbacks for UI integration
   - **Points:** 5

3. **[LOG-3] Add task timing and profiling**
   - Track per-task execution time
   - Track queue wait time
   - Generate performance reports
   - **Points:** 5

4. **[ERR-1] Enhanced error handling**
   - Aggregate errors from failed tasks
   - Add retry logic with exponential backoff
   - Partial result recovery
   - **Points:** 8

5. **[TEST-4] Error handling tests**
   - Simulate worker crashes
   - Test timeout scenarios
   - Test partial failures
   - **Points:** 5

**Total Points:** 26

---

### Sprint 5: Integration & ParallelManyBodyComputer (Weeks 9-10)
**Goal:** Complete integration with existing API

**Stories:**

1. **[INT-3] Implement ParallelManyBodyComputer.compute()**
   - Use executor to run all tasks
   - Collect results into component_results dict
   - Call existing analyze() method
   - **Points:** 8

2. **[INT-4] Add executor auto-selection**
   - Detect available parallelism
   - Choose optimal worker count
   - Environment variable configuration
   - **Points:** 3

3. **[INT-5] Update ManyBodyComputer factory**
   - Add `parallel` parameter
   - Add `n_workers` parameter
   - Add `executor_type` parameter
   - **Points:** 3

4. **[INT-6] Backward compatibility layer**
   - Ensure all existing code works unchanged
   - Add deprecation warnings for future changes
   - **Points:** 3

5. **[TEST-5] End-to-end integration tests**
   - Test full workflows with parallel execution
   - Test all BSSE types + drivers
   - Test multi-level calculations
   - **Points:** 8

**Total Points:** 25

---

### Sprint 6: Scheduler & Load Balancing (Weeks 11-12)
**Goal:** Optimize task distribution

**Stories:**

1. **[SCHED-1] Design scheduling system**
   - Priority-based scheduling
   - Cost estimation for tasks
   - Dynamic vs. static scheduling
   - **Points:** 5

2. **[SCHED-2] Implement TaskScheduler**
   - Sort tasks by priority/cost
   - Group by nbody level
   - Implement work-stealing
   - **Points:** 8

3. **[SCHED-3] Cost estimation**
   - Estimate QC calculation cost from:
     - Number of atoms
     - Method complexity (HF < MP2 < CCSD)
     - Basis set size
   - Historical timing data
   - **Points:** 5

4. **[OPT-1] Chunked submission**
   - Submit tasks in batches
   - Reduce memory footprint
   - Dynamic batch sizing
   - **Points:** 5

5. **[TEST-6] Scheduler tests**
   - Verify priority ordering
   - Test load balancing effectiveness
   - Benchmark vs. naive scheduling
   - **Points:** 5

**Total Points:** 28

---

### Sprint 7: Checkpointing & Caching (Weeks 13-14)
**Goal:** Add fault tolerance and result persistence

**Stories:**

1. **[CACHE-1] Design checkpoint format**
   - JSON or HDF5 format
   - Store task results incrementally
   - Version compatibility
   - **Points:** 3

2. **[CACHE-2] Implement checkpoint manager**
   - Save results after each batch
   - Resume from checkpoint
   - Detect completed tasks
   - **Points:** 8

3. **[CACHE-3] Result caching system**
   - Cache by task hash
   - Configurable cache location
   - TTL and eviction policies
   - **Points:** 5

4. **[CACHE-4] Integration with executors**
   - Check cache before execution
   - Save to cache after execution
   - Cache statistics
   - **Points:** 5

5. **[TEST-7] Checkpoint/cache tests**
   - Test checkpoint save/restore
   - Test cache hit/miss
   - Test corruption handling
   - **Points:** 5

**Total Points:** 26

---

### Sprint 8: MPI Foundation (Weeks 15-16)
**Goal:** Set up MPI infrastructure

**Stories:**

1. **[MPI-1] Add mpi4py dependency**
   - Make optional with extras_require
   - Add conda environment file
   - Document installation
   - **Points:** 2

2. **[MPI-2] Design MPI communication protocol**
   - Master-worker architecture
   - Message types (TASK, RESULT, SHUTDOWN)
   - Serialization format
   - **Points:** 5

3. **[MPI-3] Implement MPIExecutor skeleton**
   - Initialize MPI communicator
   - Rank-based role assignment
   - Basic send/receive
   - **Points:** 8

4. **[MPI-4] Worker loop implementation**
   - Receive tasks
   - Execute locally
   - Send results back
   - **Points:** 5

5. **[TEST-8] MPI basic tests**
   - Test with 2, 4, 8 processes
   - Local machine tests
   - **Points:** 5

**Total Points:** 25

---

### Sprint 9: MPI Advanced Features (Weeks 17-18)
**Goal:** Production-ready MPI executor

**Stories:**

1. **[MPI-5] Fault tolerance**
   - Detect worker failures
   - Reassign tasks
   - Graceful degradation
   - **Points:** 8

2. **[MPI-6] Non-blocking communication**
   - Use `Irecv`/`Isend`
   - Overlap computation and communication
   - **Points:** 5

3. **[MPI-7] HPC cluster integration**
   - SLURM job templates
   - PBS/Torque templates
   - Environment setup scripts
   - **Points:** 5

4. **[MPI-8] Performance optimization**
   - Minimize pickling overhead
   - Efficient task distribution
   - Benchmark and profile
   - **Points:** 5

5. **[TEST-9] Multi-node testing**
   - Test on real HPC cluster
   - Scaling studies (1-100 nodes)
   - Network failure simulation
   - **Points:** 8

**Total Points:** 31

---

### Sprint 10: Documentation & Release (Weeks 19-20)
**Goal:** Prepare for production release

**Stories:**

1. **[DOC-2] API documentation**
   - Docstrings for all public APIs
   - Generate Sphinx docs
   - Code examples
   - **Points:** 5

2. **[DOC-3] User tutorials**
   - Beginner: Basic parallel usage
   - Intermediate: Executor selection & tuning
   - Advanced: MPI on HPC clusters
   - **Points:** 8

3. **[DOC-4] Performance guide**
   - When to use parallel execution
   - Optimal worker count selection
   - Common pitfalls
   - **Points:** 5

4. **[DOC-5] Migration guide**
   - Converting sequential to parallel
   - API changes
   - Troubleshooting
   - **Points:** 3

5. **[REL-1] Release preparation**
   - Version bump
   - Changelog
   - PyPI package build
   - **Points:** 3

6. **[REL-2] Announcement materials**
   - Blog post
   - Example notebooks
   - Benchmark results
   - **Points:** 3

**Total Points:** 27

---

## Epic Breakdown

### Epic 1: Parallel Infrastructure 🏗️
**Duration:** Weeks 1-4
**Priority:** P0

**Description:**
Establish the foundational parallel execution framework with abstract interfaces and basic implementations.

**Components:**
- Base executor abstract class
- Task and result data models
- Sequential executor (reference implementation)
- Test infrastructure

**Success Criteria:**
- All base classes fully type-hinted
- 100% test coverage
- Documentation complete

**Dependencies:** None

---

### Epic 2: Multiprocessing Executor 🚀
**Duration:** Weeks 5-10
**Priority:** P0

**Description:**
Implement production-ready parallel execution using Python's multiprocessing module for single-node parallelism.

**Components:**
- MultiprocessingExecutor
- Task serialization
- Progress tracking and logging
- Error handling and retry logic
- Integration with ManyBodyComputer

**Success Criteria:**
- 2-4x speedup on typical workloads
- Zero deadlocks or race conditions
- Graceful error recovery

**Dependencies:** Epic 1

---

### Epic 3: Advanced Scheduling 📊
**Duration:** Weeks 11-14
**Priority:** P1

**Description:**
Optimize task distribution and add fault tolerance through intelligent scheduling and checkpointing.

**Components:**
- Priority-based task scheduler
- Cost estimation models
- Load balancing algorithms
- Checkpoint and caching system

**Success Criteria:**
- 15%+ improvement over naive scheduling
- Checkpoint/resume works reliably
- Cache hit rate >50% for repeated calculations

**Dependencies:** Epic 2

---

### Epic 4: MPI Distributed Execution 🌐
**Duration:** Weeks 15-18
**Priority:** P2

**Description:**
Enable massively parallel execution across multiple nodes using MPI for HPC environments.

**Components:**
- mpi4py integration
- Master-worker architecture
- Fault-tolerant communication
- HPC deployment tools

**Success Criteria:**
- Scales to 100+ workers
- <10% overhead vs. perfect parallelism
- Works on 3+ HPC schedulers

**Dependencies:** Epic 2, Epic 3

---

### Epic 5: Documentation & Training 📚
**Duration:** Weeks 1-20 (Continuous)
**Priority:** P1

**Description:**
Comprehensive documentation, tutorials, and training materials for users and developers.

**Components:**
- API documentation (Sphinx)
- User tutorials (Jupyter notebooks)
- Developer guide
- Performance tuning guide
- Teaching materials

**Success Criteria:**
- 95%+ docstring coverage
- Tutorials tested by external users
- Zero documentation bugs in first month

**Dependencies:** All epics (continuous)

---

## Technical Specifications

### 1. BaseParallelExecutor Interface

```python
from abc import ABC, abstractmethod
from typing import List, Optional, Callable
from dataclasses import dataclass

@dataclass
class ExecutorConfig:
    """Configuration for parallel executors."""
    n_workers: Optional[int] = None  # None = auto-detect
    timeout_per_task: float = 3600.0  # seconds
    max_retries: int = 2
    checkpoint_interval: int = 10  # tasks
    cache_dir: Optional[str] = None
    log_level: str = "INFO"

class BaseParallelExecutor(ABC):
    """Abstract base class for parallel executors.

    All executors must implement this interface to be compatible
    with ParallelManyBodyComputer.
    """

    def __init__(self, config: Optional[ExecutorConfig] = None):
        self.config = config or ExecutorConfig()
        self._is_initialized = False

    @abstractmethod
    def initialize(self) -> None:
        """Initialize executor resources (pools, communicators, etc.)."""
        pass

    @abstractmethod
    def execute(
        self,
        tasks: List[ParallelTask],
        progress_callback: Optional[Callable[[str, int, int], None]] = None
    ) -> List[TaskResult]:
        """Execute tasks in parallel.

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to execute
        progress_callback : Optional callable
            Called as progress_callback(task_id, completed, total)

        Returns
        -------
        List[TaskResult]
            Results in same order as input tasks
        """
        pass

    @abstractmethod
    def shutdown(self, wait: bool = True) -> None:
        """Shutdown executor and clean up resources."""
        pass

    def __enter__(self):
        self.initialize()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.shutdown()
        return False

    @property
    def name(self) -> str:
        """Executor name for logging."""
        return self.__class__.__name__
```

---

### 2. ParallelTask and TaskResult Models

```python
from dataclasses import dataclass, field
from typing import Optional, Dict, Any
from qcelemental.models import Molecule, AtomicInput, AtomicResult
import time

@dataclass
class ParallelTask:
    """Specification for a single parallel QC computation."""

    task_id: str  # Unique identifier (from labeler)
    chemistry: str  # Model chemistry (e.g., "mp2/cc-pvdz")
    label: str  # Full label string
    molecule: Molecule  # QCElemental molecule
    atomic_input: AtomicInput  # Full QC input specification

    # Scheduling metadata
    priority: int = 0  # Higher = execute sooner
    estimated_cost: float = 1.0  # Relative cost estimate
    nbody: int = 1  # N-body level

    # Dependencies (for future DAG scheduling)
    depends_on: List[str] = field(default_factory=list)

    def __hash__(self):
        return hash(self.task_id)

    def __lt__(self, other):
        # For priority queue
        return self.priority < other.priority

@dataclass
class TaskResult:
    """Result of a parallel task execution."""

    task_id: str
    success: bool

    # QC result (if successful)
    atomic_result: Optional[AtomicResult] = None

    # Error info (if failed)
    error_type: Optional[str] = None
    error_message: Optional[str] = None

    # Metadata
    execution_time: float = 0.0  # seconds
    attempt_number: int = 1
    worker_id: Optional[str] = None
    timestamp: float = field(default_factory=time.time)

    @property
    def return_result(self) -> Optional[float]:
        """Extract energy/gradient/hessian from result."""
        if self.success and self.atomic_result:
            return self.atomic_result.return_result
        return None
```

---

### 3. SequentialExecutor (Reference Implementation)

```python
import logging
from typing import List, Optional, Callable
from .base import BaseParallelExecutor, ParallelTask, TaskResult
from .worker import execute_single_task

logger = logging.getLogger(__name__)

class SequentialExecutor(BaseParallelExecutor):
    """Sequential executor - executes tasks one at a time.

    This executor provides no parallelism but serves as:
    1. Reference implementation for correctness testing
    2. Fallback when parallelism is unavailable
    3. Debugging tool
    """

    def initialize(self) -> None:
        """No initialization needed for sequential execution."""
        logger.info("SequentialExecutor initialized")
        self._is_initialized = True

    def execute(
        self,
        tasks: List[ParallelTask],
        progress_callback: Optional[Callable[[str, int, int], None]] = None
    ) -> List[TaskResult]:
        """Execute tasks sequentially in order."""
        results = []
        total = len(tasks)

        logger.info(f"Executing {total} tasks sequentially")

        for i, task in enumerate(tasks):
            logger.debug(f"Executing task {i+1}/{total}: {task.task_id}")

            # Execute task
            result = execute_single_task(task, self.config)
            results.append(result)

            # Progress callback
            if progress_callback:
                progress_callback(task.task_id, i+1, total)

            # Log result
            status = "✓" if result.success else "✗"
            logger.info(
                f"{status} Task {i+1}/{total} ({task.task_id}): "
                f"{result.execution_time:.2f}s"
            )

        return results

    def shutdown(self, wait: bool = True) -> None:
        """No cleanup needed."""
        logger.info("SequentialExecutor shut down")
        self._is_initialized = False
```

---

### 4. MultiprocessingExecutor (Core Implementation)

```python
import multiprocessing as mp
from multiprocessing import Pool
import logging
from typing import List, Optional, Callable
from .base import BaseParallelExecutor, ParallelTask, TaskResult, ExecutorConfig
from .worker import execute_single_task

logger = logging.getLogger(__name__)

class MultiprocessingExecutor(BaseParallelExecutor):
    """Parallel executor using Python multiprocessing.

    Uses process pool for CPU-bound QC calculations.
    Suitable for single-node parallelism.
    """

    def __init__(self, config: Optional[ExecutorConfig] = None):
        super().__init__(config)
        self._pool: Optional[Pool] = None
        self._n_workers = config.n_workers if config else None

    def initialize(self) -> None:
        """Create process pool."""
        if self._n_workers is None:
            self._n_workers = mp.cpu_count()

        logger.info(f"Initializing pool with {self._n_workers} workers")

        self._pool = Pool(
            processes=self._n_workers,
            initializer=_worker_init,
            initargs=(self.config,)
        )

        self._is_initialized = True

    def execute(
        self,
        tasks: List[ParallelTask],
        progress_callback: Optional[Callable[[str, int, int], None]] = None
    ) -> List[TaskResult]:
        """Execute tasks in parallel using process pool."""
        if not self._is_initialized:
            raise RuntimeError("Executor not initialized")

        total = len(tasks)
        logger.info(f"Executing {total} tasks on {self._n_workers} workers")

        # Submit all tasks asynchronously
        async_results = []
        for task in tasks:
            ar = self._pool.apply_async(
                execute_single_task,
                args=(task, self.config),
                error_callback=lambda e: logger.error(f"Task error: {e}")
            )
            async_results.append((task, ar))

        # Collect results with progress tracking
        results = []
        for i, (task, ar) in enumerate(async_results):
            try:
                result = ar.get(timeout=self.config.timeout_per_task)
                results.append(result)

                status = "✓" if result.success else "✗"
                logger.info(
                    f"{status} Task {i+1}/{total} ({task.task_id}): "
                    f"{result.execution_time:.2f}s"
                )

                if progress_callback:
                    progress_callback(task.task_id, i+1, total)

            except mp.TimeoutError:
                logger.error(f"Task {task.task_id} timed out")
                results.append(TaskResult(
                    task_id=task.task_id,
                    success=False,
                    error_type="TimeoutError",
                    error_message=f"Task exceeded {self.config.timeout_per_task}s"
                ))
            except Exception as e:
                logger.error(f"Task {task.task_id} failed: {e}")
                results.append(TaskResult(
                    task_id=task.task_id,
                    success=False,
                    error_type=type(e).__name__,
                    error_message=str(e)
                ))

        return results

    def shutdown(self, wait: bool = True) -> None:
        """Shutdown process pool."""
        if self._pool:
            if wait:
                self._pool.close()
                self._pool.join()
            else:
                self._pool.terminate()

            logger.info("MultiprocessingExecutor shut down")
            self._is_initialized = False

def _worker_init(config: ExecutorConfig):
    """Initialize worker process."""
    import logging
    logging.basicConfig(
        level=config.log_level,
        format='%(asctime)s - Worker %(process)d - %(levelname)s - %(message)s'
    )
```

---

### 5. ParallelManyBodyComputer Integration

```python
from typing import Optional, Dict, Any, Union
from .computer import ManyBodyComputer
from .parallel.base import BaseParallelExecutor, ExecutorConfig
from .parallel.executors import SequentialExecutor, MultiprocessingExecutor
from .parallel.task import ParallelTask, TaskResult
import logging

logger = logging.getLogger(__name__)

class ParallelManyBodyComputer(ManyBodyComputer):
    """ManyBodyComputer with parallel execution support.

    Examples
    --------
    >>> # Use default multiprocessing
    >>> computer = ParallelManyBodyComputer.from_manybodyinput(
    ...     mbin,
    ...     parallel=True,
    ...     n_workers=4
    ... )

    >>> # Use custom executor
    >>> from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig
    >>> config = ExecutorConfig(n_workers=8, timeout_per_task=7200)
    >>> executor = MultiprocessingExecutor(config)
    >>> computer = ParallelManyBodyComputer.from_manybodyinput(
    ...     mbin,
    ...     executor=executor
    ... )
    """

    def __init__(
        self,
        *args,
        executor: Optional[BaseParallelExecutor] = None,
        parallel: bool = False,
        n_workers: Optional[int] = None,
        **kwargs
    ):
        super().__init__(*args, **kwargs)

        # Set up executor
        if executor is not None:
            self.executor = executor
        elif parallel:
            config = ExecutorConfig(n_workers=n_workers)
            self.executor = MultiprocessingExecutor(config)
        else:
            self.executor = SequentialExecutor()

    @classmethod
    def from_manybodyinput(
        cls,
        input_model,
        build_tasks: bool = True,
        executor: Optional[BaseParallelExecutor] = None,
        parallel: bool = False,
        n_workers: Optional[int] = None,
    ):
        """Create computer and execute calculation with parallel support."""

        # Create computer instance (same as parent)
        computer_model = cls(
            input_data=input_model,
            executor=executor,
            parallel=parallel,
            n_workers=n_workers
        )

        if not build_tasks:
            return computer_model

        # Build task list (same as parent)
        computer_model.task_list = computer_model.build_tasks()

        # NEW: Parallel execution path
        component_properties = computer_model._execute_parallel()

        # Analysis (same as parent)
        nbody_results = computer_model.qcmb_core.analyze(
            component_properties,
            component_results=computer_model.task_list
        )

        # Format results (same as parent)
        ret = computer_model.get_results(
            nbody_results,
            computer_model.task_list
        )

        return ret

    def _execute_parallel(self) -> Dict[str, Dict[str, float]]:
        """Execute all tasks using parallel executor."""

        # Convert iterate_molecules to list of ParallelTask objects
        tasks = self._prepare_tasks()

        logger.info(f"Prepared {len(tasks)} tasks for parallel execution")

        # Execute with executor
        with self.executor as ex:
            results = ex.execute(
                tasks,
                progress_callback=self._progress_callback
            )

        # Convert TaskResult back to component_properties dict
        component_properties = self._process_results(results)

        return component_properties

    def _prepare_tasks(self) -> List[ParallelTask]:
        """Convert iterate_molecules generator to ParallelTask list."""
        tasks = []

        for chem, label, imol in self.qcmb_core.iterate_molecules():
            # Get specification for this chemistry
            spec = self.input_data.specification.specification[chem]

            # Build AtomicInput
            atomic_input = self._build_atomic_input(chem, label, imol, spec)

            # Estimate cost (simple heuristic)
            estimated_cost = self._estimate_cost(imol, chem)

            # Extract nbody from label
            _, frag, _ = delabeler(label)
            nbody = len(frag)

            task = ParallelTask(
                task_id=label,
                chemistry=chem,
                label=label,
                molecule=imol,
                atomic_input=atomic_input,
                priority=-nbody,  # Higher n-body = lower priority (do small first)
                estimated_cost=estimated_cost,
                nbody=nbody
            )
            tasks.append(task)

        return tasks

    def _process_results(self, results: List[TaskResult]) -> Dict:
        """Convert TaskResult list to component_properties dict."""
        component_properties = {}

        for result in results:
            if not result.success:
                raise RuntimeError(
                    f"Task {result.task_id} failed: {result.error_message}"
                )

            # Extract properties from AtomicResult
            ar = result.atomic_result
            props = {}

            if hasattr(ar, 'return_result') and ar.return_result is not None:
                if self.driver == "energy":
                    props["energy"] = ar.return_result
                elif self.driver == "gradient":
                    props["gradient"] = ar.return_result
                elif self.driver == "hessian":
                    props["hessian"] = ar.return_result

            # Add other properties
            if hasattr(ar, 'properties'):
                for key, val in ar.properties.dict().items():
                    if val is not None:
                        props[key] = val

            component_properties[result.task_id] = props

        return component_properties

    def _progress_callback(self, task_id: str, completed: int, total: int):
        """Log progress updates."""
        pct = 100 * completed / total
        logger.info(f"Progress: {completed}/{total} ({pct:.1f}%)")

    def _estimate_cost(self, molecule: Molecule, chemistry: str) -> float:
        """Estimate relative computational cost of task."""
        natoms = len(molecule.symbols)

        # Method scaling (very rough)
        method_costs = {
            "hf": 1.0,
            "mp2": 5.0,
            "ccsd": 20.0,
            "ccsd(t)": 100.0,
        }

        method = chemistry.split("/")[0].lower()
        method_cost = method_costs.get(method, 10.0)

        # Atom scaling (rough N^3 to N^4)
        atom_cost = natoms ** 3.5

        return method_cost * atom_cost
```

---

## Teaching & Documentation

### Tutorial 1: Introduction to Parallel Execution (Beginner)

**File:** `docs/tutorials/01_parallel_basics.md`

**Topics:**
1. Why parallelize MBE calculations?
2. Understanding task independence
3. Basic usage with `parallel=True`
4. Choosing worker count
5. Monitoring progress
6. Common pitfalls

**Code Example:**
```python
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Sequential (original)
result = ManyBodyComputer.from_manybodyinput(mbin)

# Parallel (new)
result = ManyBodyComputer.from_manybodyinput(
    mbin,
    parallel=True,
    n_workers=4  # Use 4 CPU cores
)
```

---

### Tutorial 2: Advanced Executor Configuration (Intermediate)

**File:** `docs/tutorials/02_executor_config.md`

**Topics:**
1. ExecutorConfig options
2. Custom timeout values
3. Checkpointing for long calculations
4. Result caching
5. Performance profiling
6. Troubleshooting failures

**Code Example:**
```python
from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

config = ExecutorConfig(
    n_workers=8,
    timeout_per_task=3600,  # 1 hour per task
    max_retries=3,
    checkpoint_interval=20,
    cache_dir="/scratch/mbecache"
)

executor = MultiprocessingExecutor(config)

result = ManyBodyComputer.from_manybodyinput(
    mbin,
    executor=executor
)
```

---

### Tutorial 3: HPC and MPI Execution (Advanced)

**File:** `docs/tutorials/03_mpi_clusters.md`

**Topics:**
1. MPI fundamentals for MBE
2. Installing mpi4py
3. SLURM job submission
4. Multi-node scaling
5. Debugging MPI issues
6. Performance optimization

**Code Example:**
```python
# myscript.py
from qcmanybody.parallel import MPIExecutor
from qcmanybody import ManyBodyComputer

executor = MPIExecutor()  # Auto-detects MPI environment

result = ManyBodyComputer.from_manybodyinput(
    mbin,
    executor=executor
)

# Only rank 0 writes output
if executor.rank == 0:
    print(result)
```

**SLURM Script:**
```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=4:00:00

module load python/3.10
module load openmpi/4.1

mpirun -np 64 python myscript.py
```

---

### Tutorial 4: Writing Custom Executors

**File:** `docs/tutorials/04_custom_executors.md`

**Topics:**
1. BaseParallelExecutor interface
2. Task lifecycle
3. Error handling patterns
4. Integration with external schedulers (Dask, Ray)
5. Testing custom executors

---

### Developer Guide

**File:** `docs/developer/parallel_architecture.md`

**Sections:**
1. Architecture overview
2. Adding new executor types
3. Task serialization details
4. Testing guidelines
5. Performance benchmarking
6. Contributing guidelines

---

### Performance Tuning Guide

**File:** `docs/performance/tuning.md`

**Sections:**
1. Profiling parallel execution
2. Worker count selection
   - CPU-bound: n_workers = n_cores
   - Mixed: n_workers = n_cores + 4
3. Memory considerations
4. I/O bottlenecks
5. Task granularity
6. Load balancing effectiveness
7. Case studies

**Key Recommendations:**
- Small systems (<10 fragments): May not benefit from parallelism
- Medium systems (10-50 fragments): 2-8 workers optimal
- Large systems (50+ fragments): 8-32 workers, consider MPI
- CCSD(T) calculations: Limit workers due to memory

---

## Testing Strategy

### Unit Tests

**Coverage Target:** 95%+

**Test Files:**
```
qcmanybody/parallel/tests/
├── test_base.py              # Base class interface
├── test_task_models.py       # ParallelTask, TaskResult
├── test_sequential.py        # SequentialExecutor
├── test_multiprocessing.py   # MultiprocessingExecutor
├── test_scheduler.py         # Task scheduling
├── test_checkpoint.py        # Checkpointing
├── test_mpi.py              # MPIExecutor (requires mpi4py)
└── test_integration.py       # End-to-end
```

**Key Test Scenarios:**
1. **Correctness:** Results match sequential exactly
2. **Error handling:** Worker crashes, timeouts, bad input
3. **Edge cases:** Empty task list, single task, very large molecule
4. **Concurrency:** Race conditions, deadlocks
5. **Resource cleanup:** No leaked processes/files

---

### Integration Tests

**Test against existing test suite:**
```bash
pytest qcmanybody/tests/ --executor=multiprocessing --n-workers=4
pytest qcmanybody/tests/ --executor=sequential
```

All existing tests must pass with parallel executors.

---

### Performance Tests

**Benchmarks:**

```python
# benchmarks/bench_parallel.py
import time
from qcmanybody import ManyBodyComputer

def benchmark_executor(mbin, executor_type, n_workers):
    start = time.time()
    result = ManyBodyComputer.from_manybodyinput(
        mbin,
        parallel=(executor_type != "sequential"),
        n_workers=n_workers
    )
    elapsed = time.time() - start
    return elapsed

# Test systems:
# - Small: He3 trimer, 3 fragments
# - Medium: (H2O)5, 5 fragments
# - Large: (H2O)10, 10 fragments

# Configurations:
# - Sequential
# - Multiprocessing: 2, 4, 8 workers
# - MPI: 4, 8, 16, 32 workers
```

**Acceptance Criteria:**
- 2x speedup with 4 workers (medium system)
- 3x speedup with 8 workers (large system)
- <20% overhead vs. theoretical maximum

---

### Regression Tests

**Continuous Integration:**
```yaml
# .github/workflows/parallel_tests.yml
name: Parallel Module Tests

on: [push, pull_request]

jobs:
  test-parallel:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, 3.10, 3.11]
        executor: [sequential, multiprocessing]

    steps:
      - uses: actions/checkout@v3
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}
      - name: Install dependencies
        run: |
          pip install -e .[parallel]
          pip install pytest pytest-cov
      - name: Run tests
        run: |
          pytest qcmanybody/parallel/tests/ \
            --executor=${{ matrix.executor }} \
            --cov=qcmanybody.parallel \
            --cov-report=xml
      - name: Upload coverage
        uses: codecov/codecov-action@v3
```

---

## Risk Management

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| **Serialization failures** | Medium | High | Extensive testing of molecule/input pickling; fallback to JSON |
| **Deadlocks in multiprocessing** | Low | High | Use timeouts; automated deadlock detection tests |
| **MPI installation complexity** | High | Medium | Make MPI optional; provide conda environments |
| **Memory exhaustion with many workers** | Medium | High | Add memory monitoring; auto-reduce worker count |
| **QCEngine thread safety** | Low | Critical | Verify with QCEngine maintainers; add locks if needed |
| **Performance regression** | Low | Medium | Automated benchmarks in CI |

### Project Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| **Scope creep** | Medium | Medium | Strict milestone gates; defer non-critical features |
| **Resource availability** | Low | High | Prioritize Epics 1-2; defer MPI if needed |
| **Breaking API changes** | Low | High | Maintain backward compatibility; extensive deprecation period |
| **Adoption resistance** | Medium | Low | Focus on documentation; provide migration guide |

---

## Success Metrics

### Performance Metrics

1. **Speedup Ratio**
   - Target: 0.8 * n_workers (80% efficiency)
   - Measurement: Total time (sequential) / Total time (parallel with N workers)

2. **Task Throughput**
   - Target: >90% worker utilization
   - Measurement: (Worker active time) / (Total wall time * n_workers)

3. **Overhead**
   - Target: <5% for systems with >20 tasks
   - Measurement: (Parallel time - Sum of task times) / Parallel time

### Quality Metrics

1. **Code Coverage**
   - Target: >95% for parallel module
   - Tool: pytest-cov

2. **Documentation Coverage**
   - Target: 100% public APIs have docstrings
   - Tool: interrogate

3. **Bug Density**
   - Target: <0.1 bugs per 1000 lines of code
   - Measurement: GitHub issues labeled "bug"

### Adoption Metrics

1. **Usage**
   - Target: 50% of users try parallel execution in first quarter
   - Measurement: Telemetry (opt-in)

2. **User Satisfaction**
   - Target: >4.0/5.0 in user survey
   - Measurement: Post-release survey

---

## Appendix A: Task Dependencies

### Sprint Dependency Graph

```
Sprint 1 (Setup)
    ↓
Sprint 2 (Sequential)
    ↓
Sprint 3 (Multiprocessing Core) ← Sprint 4 (Logging)
    ↓                                      ↓
Sprint 5 (Integration) ←-------------------+
    ↓
Sprint 6 (Scheduling)
    ↓
Sprint 7 (Checkpointing)
    ↓
Sprint 8 (MPI Foundation)
    ↓
Sprint 9 (MPI Advanced)
    ↓
Sprint 10 (Documentation)
```

---

## Appendix B: Code Style Guide

**Standards:**
- PEP 8 compliant
- Type hints required for all public APIs
- Google-style docstrings
- Maximum line length: 100 characters
- Black formatter
- isort for imports

**Example:**
```python
from typing import List, Optional
import numpy as np

def process_results(
    results: List[TaskResult],
    tolerance: float = 1e-8
) -> Optional[np.ndarray]:
    """Process task results and compute final array.

    Parameters
    ----------
    results : List[TaskResult]
        Results from parallel execution
    tolerance : float, optional
        Numerical tolerance for validation, by default 1e-8

    Returns
    -------
    Optional[np.ndarray]
        Processed array, or None if validation fails

    Raises
    ------
    ValueError
        If any results are invalid
    """
    pass
```

---

## Appendix C: Glossary

- **BSSE:** Basis Set Superposition Error
- **MBE:** Many-Body Expansion
- **CP:** Counterpoise correction
- **NOCP:** No counterpoise
- **VMFC:** Valiron-Mayer Function Counterpoise
- **HPC:** High-Performance Computing
- **MPI:** Message Passing Interface
- **QC:** Quantum Chemistry
- **Executor:** Object responsible for parallel task execution
- **Worker:** Process/thread executing tasks
- **Task:** Single QC calculation to be executed
- **Subsystem:** Subset of fragments in the full system

---

## Appendix D: References

1. QCEngine Documentation: https://github.com/MolSSI/QCEngine
2. QCElemental Documentation: https://github.com/MolSSI/QCElemental
3. Python multiprocessing: https://docs.python.org/3/library/multiprocessing.html
4. mpi4py Documentation: https://mpi4py.readthedocs.io/
5. Many-Body Expansion Theory: J. Chem. Phys. 137, 164104 (2012)

---

## Current Status Summary (2025-11-11)

### ✅ Completed Work
1. **Parallel Module Infrastructure**
   - Base executor interface (`BaseParallelExecutor`, `ExecutorConfig`)
   - Task models (`ParallelTask`, `TaskResult`)
   - Sequential executor (reference implementation)
   - Multiprocessing executor (core functionality)
   - ParallelManyBodyComputer class
   - Worker execution functions
   - Utility functions

2. **CLI Merge**
   - Complete CLI system merged into main
   - Input file parsing (JSON/YAML)
   - Commands: run, validate, plan, convert
   - Input schema for molecule, calculation, BSSE, output

### 🎯 Current Sprint: CLI Integration
**Priority Tasks:**
1. Add `ExecutionSchema` to CLI input schema
2. Update `run.py` to use `ParallelManyBodyComputer`
3. Add `--parallel` and `--n-workers` CLI flags
4. Create parallel execution examples
5. Write integration tests
6. Update documentation

### 📋 Next Steps (In Order)
1. **Implement CLI Integration (Sprint 2.5)**
   - Add execution configuration to input schema
   - Modify run command to support parallel execution
   - Add CLI arguments for parallel options
   - Create examples and documentation

2. **Complete Testing & Validation (Milestone 4)**
   - Unit tests for all executors
   - Integration tests with CLI
   - Performance benchmarking
   - Memory leak testing

3. **Advanced Features (Milestone 5)**
   - Load balancing
   - Checkpointing
   - Result caching
   - Performance optimization

4. **MPI Support (Milestone 6)** - Future
5. **Production Release (Milestone 7)** - Future

### 🔗 Dependencies & Blockers
**No blockers.** CLI and parallel execution have been cleanly merged. Ready to proceed with integration.

### 📊 Progress Metrics
- **Milestones Complete:** 2 of 7 (29%)
- **Sprints Complete:** 2 of 10 (20%)
- **Core Infrastructure:** 90% complete
- **CLI Integration:** 0% complete (next focus)
- **Testing Coverage:** ~10% (needs expansion)

---

**Document Version:** 1.1
**Last Updated:** 2025-11-11
**Next Review:** End of Sprint 2.5 (CLI Integration)
**Maintained By:** QCManyBody Development Team
