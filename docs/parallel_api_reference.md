# Parallel Execution API Reference

Complete API documentation for QCManyBody parallel execution module.

## Module Overview

```python
from qcmanybody.parallel import (
    # Core classes
    ParallelManyBodyComputer,
    BaseParallelExecutor,
    ExecutorConfig,

    # Executors
    SequentialExecutor,
    MultiprocessingExecutor,

    # Task models
    ParallelTask,
    TaskResult,
    TaskStatus,

    # Utilities
    estimate_task_cost,
    assign_priorities,
    compute_execution_stats,
    parallel_compute_from_manybodyinput,
)
```

---

## ParallelManyBodyComputer

Main class for parallel many-body calculations.

### Class Methods

#### `from_manybodyinput()`

Create and execute parallel many-body calculation from input.

```python
@classmethod
def from_manybodyinput(
    cls,
    input_model: ManyBodyInput,
    parallel: bool = False,
    n_workers: Optional[int] = None,
    executor: Optional[BaseParallelExecutor] = None,
    build_tasks: bool = True,
    return_metadata: bool = False
) -> Union[ManyBodyResult, Dict]
```

**Parameters:**

- **input_model** (*ManyBodyInput*): Input specification for calculation
- **parallel** (*bool*, optional): Enable parallel execution. Default: `False`
- **n_workers** (*int*, optional): Number of worker processes. `None` = auto-detect. Default: `None`
- **executor** (*BaseParallelExecutor*, optional): Custom executor instance. If provided, overrides `parallel` and `n_workers`. Default: `None`
- **build_tasks** (*bool*, optional): Whether to build and execute tasks. Default: `True`
- **return_metadata** (*bool*, optional): Return execution statistics. Default: `False`

**Returns:**

- **ManyBodyResult**: Calculation results (if `return_metadata=False`)
- **Dict**: `{"result": ManyBodyResult, "metadata": Dict}` (if `return_metadata=True`)

**Raises:**

- **ValueError**: Invalid input configuration
- **RuntimeError**: Execution failure
- **ImportError**: Missing dependencies (qcengine)

**Examples:**

```python
# Simple parallel execution
result = ParallelManyBodyComputer.from_manybodyinput(
    input_model,
    parallel=True,
    n_workers=4
)

# With custom executor
from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

config = ExecutorConfig(n_workers=8, timeout_per_task=3600.0)
executor = MultiprocessingExecutor(config)

result = ParallelManyBodyComputer.from_manybodyinput(
    input_model,
    executor=executor
)

# With metadata
output = ParallelManyBodyComputer.from_manybodyinput(
    input_model,
    parallel=True,
    return_metadata=True
)

print(f"Total time: {output['metadata']['total_time']:.2f}s")
print(f"Speedup: {output['metadata']['speedup']:.2f}x")
```

---

## ExecutorConfig

Configuration for task executors.

### Constructor

```python
@dataclass
class ExecutorConfig:
    n_workers: Optional[int] = None
    timeout_per_task: float = 3600.0
    max_retries: int = 2
    checkpoint_interval: int = 10
    checkpoint_file: Optional[str] = None
    cache_dir: Optional[str] = None
    log_level: str = "INFO"
```

**Parameters:**

- **n_workers** (*int*, optional): Number of worker processes
  - `None` = auto-detect CPU count
  - Must be ≥ 1 if specified
  - Default: `None`

- **timeout_per_task** (*float*): Maximum time per task in seconds
  - Must be > 0
  - Tasks exceeding timeout are terminated
  - Default: `3600.0` (1 hour)

- **max_retries** (*int*): Maximum retry attempts for failed tasks
  - Must be ≥ 0
  - `0` = no retries, fail immediately
  - Default: `2`

- **checkpoint_interval** (*int*): Save checkpoint every N tasks
  - Must be ≥ 1
  - Currently unused (future feature)
  - Default: `10`

- **checkpoint_file** (*str*, optional): Path to checkpoint file
  - Currently unused (future feature)
  - Default: `None`

- **cache_dir** (*str*, optional): Directory for result caching
  - Currently unused (future feature)
  - Default: `None`

- **log_level** (*str*): Logging level
  - Valid: "DEBUG", "INFO", "WARNING", "ERROR"
  - Default: `"INFO"`

**Validation:**

All parameters are validated on construction. Invalid values raise `ValueError`.

**Examples:**

```python
# Default configuration
config = ExecutorConfig()

# Production configuration
config = ExecutorConfig(
    n_workers=8,
    timeout_per_task=7200.0,  # 2 hours
    max_retries=3,
    log_level="INFO"
)

# Debug configuration
config = ExecutorConfig(
    n_workers=1,
    timeout_per_task=60.0,
    max_retries=0,  # Fail fast
    log_level="DEBUG"
)
```

---

## BaseParallelExecutor

Abstract base class for all executors.

### Abstract Methods

Subclasses must implement:

```python
def initialize(self) -> None:
    """Initialize executor resources."""

def execute(
    self,
    tasks: List[ParallelTask],
    progress_callback: Optional[Callable] = None
) -> List[TaskResult]:
    """Execute tasks and return results."""

def shutdown(self, wait: bool = True) -> None:
    """Clean up executor resources."""
```

### Properties

```python
@property
def name(self) -> str:
    """Executor name (e.g., 'SequentialExecutor')."""

@property
def is_initialized(self) -> bool:
    """Whether executor is initialized and ready."""
```

### Methods

#### `validate_tasks()`

Validate task list before execution.

```python
def validate_tasks(self, tasks: List[ParallelTask]) -> None
```

**Raises:**
- **ValueError**: Empty task list
- **ValueError**: Duplicate task IDs

#### `get_info()`

Get executor information.

```python
def get_info(self) -> Dict[str, Any]
```

**Returns:**

```python
{
    "name": "SequentialExecutor",
    "n_workers": 1,
    "is_initialized": True,
    "config": ExecutorConfig(...)
}
```

### Context Manager

All executors support context manager protocol:

```python
with executor:
    results = executor.execute(tasks)
# Automatically calls initialize() and shutdown()
```

---

## SequentialExecutor

Reference implementation that executes tasks sequentially.

### Constructor

```python
def __init__(self, config: Optional[ExecutorConfig] = None)
```

**Parameters:**

- **config** (*ExecutorConfig*, optional): Configuration. If `None`, uses defaults. `n_workers` is ignored (always 1).

### Usage

```python
from qcmanybody.parallel import SequentialExecutor, ExecutorConfig

config = ExecutorConfig(
    timeout_per_task=600.0,
    log_level="DEBUG"
)

executor = SequentialExecutor(config)

with executor:
    results = executor.execute(tasks)

# Verify results
for result in results:
    if result.success:
        print(f"{result.task_id}: {result.return_result}")
    else:
        print(f"{result.task_id}: FAILED - {result.error_message}")
```

**Characteristics:**

- ✅ Simple, reliable, deterministic
- ✅ Easy to debug
- ✅ Low memory overhead
- ❌ No parallelization benefit
- ❌ Slow for large calculations

**Use cases:**

- Development and debugging
- Small calculations (< 10 tasks)
- Systems with single CPU core
- Baseline for performance comparison

---

## MultiprocessingExecutor

Production executor using Python multiprocessing for parallel execution.

### Constructor

```python
def __init__(self, config: Optional[ExecutorConfig] = None)
```

**Parameters:**

- **config** (*ExecutorConfig*, optional): Configuration. If `None`, uses defaults with auto-detected CPU count.

### Usage

```python
from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

config = ExecutorConfig(
    n_workers=4,
    timeout_per_task=1800.0,
    max_retries=2
)

executor = MultiprocessingExecutor(config)

with executor:
    results = executor.execute(tasks)
```

**Implementation Details:**

- Uses `multiprocessing.Pool` for worker management
- Each worker is a separate process (true parallelism)
- Tasks are submitted asynchronously to pool
- Results are collected in original task order
- Workers are initialized once and reused

**Characteristics:**

- ✅ True parallelism (multi-core)
- ✅ Scales with CPU cores
- ✅ Isolated worker processes (crash-safe)
- ❌ Process creation overhead
- ❌ Requires pickleable data

**Use cases:**

- Production calculations
- Multi-core systems
- Medium to large calculations (10-1000 tasks)

### Special Behaviors

**Auto-detection:**

If `config.n_workers` is `None`, detects available CPU cores:

```python
import multiprocessing as mp
n_workers = mp.cpu_count()  # Physical + logical cores
```

**Error handling:**

- Failed tasks return `TaskResult` with `success=False`
- Timeout tasks are terminated and marked as `TIMEOUT`
- Worker crashes are caught and reported

**Resource cleanup:**

- `shutdown()` waits for workers to finish by default
- `shutdown(wait=False)` terminates workers immediately
- Pool is always cleaned up, even on exception

---

## ParallelTask

Specification for a single quantum chemistry calculation task.

### Constructor

```python
@dataclass
class ParallelTask:
    task_id: str
    chemistry: str
    label: str
    molecule: Molecule
    atomic_input: AtomicInput
    priority: int = 0
    estimated_cost: float = 1.0
    nbody: int = 1
    depends_on: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
```

**Parameters:**

- **task_id** (*str*): Unique identifier for task
- **chemistry** (*str*): Chemistry specification (e.g., "hf/sto-3g")
- **label** (*str*): Human-readable label
- **molecule** (*Molecule*): Molecular geometry
- **atomic_input** (*AtomicInput*): QCEngine input specification
- **priority** (*int*, optional): Task priority (higher = more important). Default: `0`
- **estimated_cost** (*float*, optional): Estimated computational cost. Default: `1.0`
- **nbody** (*int*, optional): N-body level (1, 2, 3, ...). Default: `1`
- **depends_on** (*List[str]*, optional): Task IDs of dependencies. Default: `[]`
- **metadata** (*Dict*, optional): Additional metadata. Default: `{}`

### Methods

#### Comparison

Tasks support comparison for priority queue ordering:

```python
task1 = ParallelTask(..., priority=10)
task2 = ParallelTask(..., priority=5)

assert task1 < task2  # Higher priority comes first
```

#### Hash and Equality

Tasks are hashable by `task_id`:

```python
task1 = ParallelTask(task_id="task_1", ...)
task2 = ParallelTask(task_id="task_1", ...)  # Same ID

assert task1 == task2
assert hash(task1) == hash(task2)

task_set = {task1, task2}
assert len(task_set) == 1  # Deduplicated
```

---

## TaskResult

Result of task execution.

### Constructor

```python
@dataclass
class TaskResult:
    task_id: str
    success: bool
    status: TaskStatus = TaskStatus.PENDING
    atomic_result: Optional[AtomicResult] = None
    error_type: Optional[str] = None
    error_message: Optional[str] = None
    error_traceback: Optional[str] = None
    execution_time: float = 0.0
    queue_time: float = 0.0
    worker_id: Optional[str] = None
    attempt_number: int = 1
    metadata: Dict[str, Any] = field(default_factory=dict)
```

**Parameters:**

- **task_id** (*str*): Task identifier (matches ParallelTask.task_id)
- **success** (*bool*): Whether task succeeded
- **status** (*TaskStatus*): Task status (COMPLETED, FAILED, TIMEOUT, etc.)
- **atomic_result** (*AtomicResult*, optional): QCEngine result if successful
- **error_type** (*str*, optional): Exception class name if failed
- **error_message** (*str*, optional): Error message if failed
- **error_traceback** (*str*, optional): Full traceback if failed
- **execution_time** (*float*): Time spent executing (seconds)
- **queue_time** (*float*): Time spent in queue (seconds)
- **worker_id** (*str*, optional): Worker that executed task
- **attempt_number** (*int*): Attempt number (1 = first try)
- **metadata** (*Dict*, optional): Additional metadata

### Properties

```python
@property
def return_result(self) -> Optional[float]:
    """Extract return_result from atomic_result, or None if failed."""

@property
def total_time(self) -> float:
    """Total time = execution_time + queue_time."""
```

### Methods

```python
def to_dict(self) -> Dict[str, Any]:
    """Convert to dictionary for serialization."""
```

### String Representation

```python
result = TaskResult(task_id="task_1", success=True, status=TaskStatus.COMPLETED, execution_time=5.2)
print(result)
# Output: TaskResult(✓ task_id='task_1', status=completed, time=5.20s, attempt=1)
```

---

## TaskStatus

Enumeration of task statuses.

```python
class TaskStatus(Enum):
    PENDING = "pending"        # Not yet started
    RUNNING = "running"        # Currently executing
    COMPLETED = "completed"    # Finished successfully
    FAILED = "failed"          # Failed with error
    TIMEOUT = "timeout"        # Exceeded timeout
    CANCELLED = "cancelled"    # Cancelled by user
```

---

## Utility Functions

### estimate_task_cost()

Estimate computational cost of a task.

```python
def estimate_task_cost(task: ParallelTask) -> float
```

**Returns:** Estimated relative cost (1.0 = baseline)

**Heuristic:**

```python
cost = (
    n_atoms ** 2.5 *
    method_factor *
    basis_factor
)
```

Where:
- `method_factor`: HF=1.0, MP2=5.0, CCSD=20.0
- `basis_factor`: STO-3G=1.0, cc-pVDZ=5.0, cc-pVTZ=15.0

**Example:**

```python
cost = estimate_task_cost(task)
print(f"Estimated cost: {cost:.1f}x baseline")
```

### assign_priorities()

Assign priorities to tasks based on strategy.

```python
def assign_priorities(
    tasks: List[ParallelTask],
    strategy: str = "nbody_first"
) -> None
```

**Parameters:**

- **tasks** (*List[ParallelTask]*): Tasks to prioritize (modified in-place)
- **strategy** (*str*): Priority strategy
  - `"nbody_first"`: Higher n-body levels first (default)
  - `"cost_first"`: Higher cost tasks first
  - `"fifo"`: First-in-first-out (no prioritization)

**Example:**

```python
assign_priorities(tasks, strategy="nbody_first")

# Now tasks are prioritized:
# - 3-body tasks: priority = 3
# - 2-body tasks: priority = 2
# - 1-body tasks: priority = 1
```

### compute_execution_stats()

Compute statistics from execution results.

```python
def compute_execution_stats(results: List[TaskResult]) -> Dict[str, Any]
```

**Returns:**

```python
{
    "total_tasks": 14,
    "successful": 14,
    "failed": 0,
    "timeout": 0,
    "total_time": 45.3,
    "avg_task_time": 3.2,
    "min_task_time": 1.5,
    "max_task_time": 5.8,
    "speedup": 3.7,  # If sequential time available
}
```

**Example:**

```python
stats = compute_execution_stats(results)
print(f"Success rate: {stats['successful']/stats['total_tasks']*100:.1f}%")
print(f"Average time: {stats['avg_task_time']:.2f}s per task")
```

### parallel_compute_from_manybodyinput()

Convenience function for parallel execution.

```python
def parallel_compute_from_manybodyinput(
    input_model: ManyBodyInput,
    **kwargs
) -> ManyBodyResult
```

Equivalent to:

```python
ParallelManyBodyComputer.from_manybodyinput(input_model, **kwargs)
```

---

## Type Hints

The module provides comprehensive type hints for type checking:

```python
from typing import List, Optional, Dict, Any, Callable
from qcmanybody.parallel import (
    ParallelTask,
    TaskResult,
    BaseParallelExecutor,
)

def process_results(
    results: List[TaskResult],
    callback: Optional[Callable[[TaskResult], None]] = None
) -> Dict[str, Any]:
    """Process task results with type safety."""
    stats: Dict[str, Any] = {}

    for result in results:
        if result.success:
            if callback:
                callback(result)
        else:
            print(f"Failed: {result.error_message}")

    return stats
```

---

## Error Handling

### Common Exceptions

**ValueError:**
- Invalid configuration (e.g., `n_workers=0`)
- Empty task list
- Duplicate task IDs

**RuntimeError:**
- Executor not initialized
- Execution failure
- Double initialization

**ImportError:**
- Missing qcengine dependency
- Missing QC program

**TimeoutError:**
- Task exceeds `timeout_per_task` (caught and returned as TaskResult)

### Best Practices

```python
from qcmanybody.parallel import (
    MultiprocessingExecutor,
    ExecutorConfig,
    ParallelManyBodyComputer,
)

try:
    config = ExecutorConfig(
        n_workers=4,
        timeout_per_task=1800.0,
        max_retries=2,
        log_level="INFO"
    )

    executor = MultiprocessingExecutor(config)

    result = ParallelManyBodyComputer.from_manybodyinput(
        input_model,
        executor=executor
    )

except ValueError as e:
    print(f"Configuration error: {e}")

except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Install with: pip install qcengine")

except RuntimeError as e:
    print(f"Execution failed: {e}")
    # Check logs for details

finally:
    # Executor cleanup handled by context manager
    pass
```

---

## Thread Safety

**Executor classes** are NOT thread-safe. Do not share executor instances across threads.

**Safe:**

```python
# Each thread has its own executor
def worker():
    executor = MultiprocessingExecutor(config)
    with executor:
        results = executor.execute(tasks)

threads = [Thread(target=worker) for _ in range(2)]
```

**Unsafe:**

```python
# Shared executor across threads - DON'T DO THIS
executor = MultiprocessingExecutor(config)

def worker():
    with executor:  # Race condition!
        results = executor.execute(tasks)
```

---

## Performance Considerations

### Memory Usage

Each executor has different memory characteristics:

**SequentialExecutor:**
- Memory: 1x QC program + 1x task
- Predictable, low memory usage

**MultiprocessingExecutor:**
- Memory: N × (QC program + task)
- Where N = `n_workers`
- Can be significant for large systems

### CPU Usage

**SequentialExecutor:**
- CPU: 1 core at 100%
- Other cores idle

**MultiprocessingExecutor:**
- CPU: N cores at ~100%
- Scales well up to physical core count
- Hyperthreading provides minimal benefit

### I/O Considerations

All executors may bottleneck on:
- Disk I/O (reading/writing scratch files)
- Network I/O (accessing remote filesystems)

**Mitigation:**
- Use local scratch directories
- Reduce `n_workers` if I/O is bottleneck
- Consider SSD for scratch files

---

## See Also

- **User Guide:** `docs/parallel_execution_guide.md`
- **Development Plan:** `PARALLEL_EXECUTION_PLAN.md`
- **Module README:** `qcmanybody/parallel/README.md`
- **Examples:** `examples/cli/05_parallel_basic.json`
