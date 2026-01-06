# Parallel Execution API Design

## Overview

This document defines the API design for parallel execution features in QCManyBody. The design prioritizes backward compatibility while providing powerful new parallel execution capabilities.

## Design Principles

1. **Backward Compatibility**: Existing code continues to work unchanged
2. **Progressive Enhancement**: Users can opt-in to parallel features
3. **Configuration Flexibility**: Multiple ways to configure parallel execution
4. **Fail-Safe Defaults**: Automatic fallback to sequential execution
5. **Resource Awareness**: Automatic adaptation to system capabilities

## Core API Components

### 1. ParallelManyBodyComputer

The main entry point for parallel execution, extending the existing `ManyBodyComputer` interface.

```python
from qcmanybody import ParallelManyBodyComputer, ParallelMode

class ParallelManyBodyComputer(ManyBodyComputer):
    """Enhanced ManyBodyComputer with parallel execution capabilities."""

    def __init__(self,
                 input_data: ManyBodyInput,
                 parallel_config: Optional[ParallelConfig] = None,
                 **kwargs):
        """
        Initialize parallel many-body computer.

        Parameters
        ----------
        input_data : ManyBodyInput
            Standard QCManyBody input specification
        parallel_config : ParallelConfig, optional
            Parallel execution configuration
        **kwargs
            Additional arguments passed to ManyBodyComputer
        """
```

### 2. ParallelConfig

Configuration class for parallel execution parameters.

```python
@dataclass
class ParallelConfig:
    """Configuration for parallel execution."""

    # Execution mode
    mode: ParallelMode = ParallelMode.MULTIPROCESSING
    max_workers: Optional[int] = None  # Auto-detect if None

    # Resource management
    memory_limit_per_worker: Optional[str] = None  # e.g., "2GB"
    timeout_per_fragment: Optional[int] = None     # seconds

    # Load balancing
    load_balancing: LoadBalanceStrategy = LoadBalanceStrategy.DYNAMIC
    work_stealing: bool = True

    # Error handling
    retry_failed_fragments: bool = True
    max_retries: int = 2
    fallback_to_sequential: bool = True

    # Progress monitoring
    progress_reporting: bool = True
    log_level: str = "INFO"
```

### 3. Parallel Execution Modes

```python
class ParallelMode(Enum):
    """Available parallel execution modes."""

    SEQUENTIAL = "sequential"        # No parallelization
    MULTIPROCESSING = "multiprocessing"  # Process-based parallelization
    MPI = "mpi"                     # MPI-based distributed execution
    ADAPTIVE = "adaptive"           # Auto-select best mode
```

### 4. Load Balancing Strategies

```python
class LoadBalanceStrategy(Enum):
    """Load balancing strategies for parallel execution."""

    STATIC = "static"               # Pre-distribute work evenly
    DYNAMIC = "dynamic"             # Distribute based on cost estimation
    WORK_STEALING = "work_stealing" # Workers steal from busy workers
    ADAPTIVE = "adaptive"           # Adapt strategy based on performance
```

## Usage Examples

### Basic Parallel Execution

```python
from qcmanybody import ParallelManyBodyComputer, ParallelConfig, ParallelMode

# Simple parallel execution with default settings
computer = ParallelManyBodyComputer(
    input_data=mbe_input,
    parallel_config=ParallelConfig(mode=ParallelMode.MULTIPROCESSING)
)

result = computer.compute()
```

### Advanced Configuration

```python
# Advanced parallel configuration
parallel_config = ParallelConfig(
    mode=ParallelMode.MULTIPROCESSING,
    max_workers=8,
    memory_limit_per_worker="4GB",
    load_balancing=LoadBalanceStrategy.DYNAMIC,
    timeout_per_fragment=3600,  # 1 hour per fragment
    retry_failed_fragments=True,
    progress_reporting=True
)

computer = ParallelManyBodyComputer(
    input_data=mbe_input,
    parallel_config=parallel_config
)

result = computer.compute()
```

### MPI Execution

```python
# MPI-based execution for HPC clusters
parallel_config = ParallelConfig(
    mode=ParallelMode.MPI,
    max_workers=16,  # Will use MPI communicator size if not specified
    load_balancing=LoadBalanceStrategy.STATIC,
    memory_limit_per_worker="2GB"
)

computer = ParallelManyBodyComputer(
    input_data=mbe_input,
    parallel_config=parallel_config
)

result = computer.compute()
```

### Environment-Based Configuration

```python
# Configuration via environment variables
import os
os.environ['QCMB_PARALLEL_MODE'] = 'multiprocessing'
os.environ['QCMB_MAX_WORKERS'] = '4'
os.environ['QCMB_MEMORY_LIMIT'] = '2GB'

# Auto-configure from environment
computer = ParallelManyBodyComputer.from_environment(input_data=mbe_input)
result = computer.compute()
```

## Backward Compatibility

### Existing Code Unchanged

```python
# This continues to work exactly as before
from qcmanybody import ManyBodyComputer

computer = ManyBodyComputer(input_data=mbe_input)
result = computer.compute()  # Sequential execution as always
```

### Progressive Enhancement

```python
# Easy migration path - just change the class
from qcmanybody import ParallelManyBodyComputer

# No other changes needed - defaults to sequential execution
computer = ParallelManyBodyComputer(input_data=mbe_input)
result = computer.compute()

# Enable parallel execution with minimal changes
computer = ParallelManyBodyComputer(
    input_data=mbe_input,
    parallel_config={'mode': 'multiprocessing'}
)
result = computer.compute()
```

## Error Handling and Fallbacks

### Automatic Fallback

```python
class ParallelManyBodyComputer:
    def compute(self) -> ManyBodyResult:
        """Compute with automatic fallback to sequential execution."""
        try:
            return self._compute_parallel()
        except ParallelExecutionError as e:
            if self.parallel_config.fallback_to_sequential:
                logger.warning(f"Parallel execution failed: {e}")
                logger.info("Falling back to sequential execution")
                return self._compute_sequential()
            else:
                raise
```

### Retry Logic

```python
class FragmentExecutor:
    def execute_fragment(self, fragment) -> AtomicResult:
        """Execute fragment with retry logic."""
        last_exception = None

        for attempt in range(self.max_retries + 1):
            try:
                return self._execute_single(fragment)
            except Exception as e:
                last_exception = e
                if attempt < self.max_retries:
                    logger.warning(f"Fragment failed (attempt {attempt + 1}): {e}")
                    continue

        # All retries exhausted
        raise FragmentExecutionError(f"Fragment failed after {self.max_retries} retries") from last_exception
```

## Resource Management

### Memory Management

```python
class ResourceManager:
    def __init__(self, config: ParallelConfig):
        self.config = config
        self.available_memory = self._get_available_memory()

    def get_optimal_worker_count(self) -> int:
        """Determine optimal worker count based on available resources."""
        if self.config.max_workers:
            requested_workers = self.config.max_workers
        else:
            requested_workers = os.cpu_count()

        if self.config.memory_limit_per_worker:
            memory_constrained_workers = self.available_memory // self._parse_memory_limit(self.config.memory_limit_per_worker)
            return min(requested_workers, memory_constrained_workers)

        return requested_workers
```

### Process Lifecycle Management

```python
class WorkerPool:
    def __enter__(self):
        """Initialize worker processes."""
        self.workers = self._create_workers()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Clean up worker processes."""
        self._shutdown_workers(graceful=True)
        self._cleanup_resources()
```

## Progress Monitoring

### Progress Reporting Interface

```python
class ProgressReporter:
    def on_level_started(self, level: int, fragment_count: int):
        """Called when starting a new N-body level."""

    def on_fragment_completed(self, fragment_label: str, success: bool):
        """Called when a fragment calculation completes."""

    def on_level_completed(self, level: int, duration: float):
        """Called when an N-body level completes."""

    def on_calculation_completed(self, total_duration: float):
        """Called when entire calculation completes."""
```

### Usage with Progress Reporting

```python
from qcmanybody.progress import ConsoleProgressReporter, FileProgressReporter

# Console progress reporting
progress_reporter = ConsoleProgressReporter(show_progress_bar=True)

computer = ParallelManyBodyComputer(
    input_data=mbe_input,
    parallel_config=ParallelConfig(
        mode=ParallelMode.MULTIPROCESSING,
        progress_reporting=True
    ),
    progress_reporter=progress_reporter
)

result = computer.compute()
```

## Integration Points

### QCEngine Integration

```python
class QCEngineAdapter:
    """Adapter for QCEngine parallel execution."""

    def __init__(self, qc_program: str, config: ParallelConfig):
        self.qc_program = qc_program
        self.config = config

    def execute_fragment(self, atomic_input: AtomicInput) -> AtomicResult:
        """Execute single fragment via QCEngine."""
        return qcng.compute(
            atomic_input,
            self.qc_program,
            raise_error=False
        )
```

### Task Distribution

```python
class TaskDistributor:
    """Distribute fragments across workers."""

    def distribute_level(self,
                        fragments: List[FragmentTask],
                        workers: int) -> List[List[FragmentTask]]:
        """Distribute fragments for a given N-body level."""

        if self.load_balance_strategy == LoadBalanceStrategy.STATIC:
            return self._static_distribution(fragments, workers)
        elif self.load_balance_strategy == LoadBalanceStrategy.DYNAMIC:
            return self._dynamic_distribution(fragments, workers)
        else:
            raise ValueError(f"Unknown load balance strategy: {self.load_balance_strategy}")
```

## Testing and Validation

### Test Utilities

```python
class ParallelTestSuite:
    """Utilities for testing parallel execution."""

    def validate_parallel_results(self,
                                parallel_result: ManyBodyResult,
                                sequential_result: ManyBodyResult,
                                tolerance: float = 1e-10) -> bool:
        """Validate that parallel and sequential results are identical."""

    def benchmark_performance(self,
                            system: TestSystem,
                            parallel_config: ParallelConfig) -> BenchmarkResult:
        """Benchmark parallel execution performance."""
```

## Future Extensions

### Plugin Architecture

```python
class ParallelExecutionPlugin:
    """Base class for parallel execution plugins."""

    def can_handle(self, config: ParallelConfig) -> bool:
        """Check if this plugin can handle the given configuration."""

    def create_executor(self, config: ParallelConfig) -> ParallelExecutor:
        """Create parallel executor for this configuration."""
```

### Advanced Load Balancing

```python
class MLLoadBalancer:
    """Machine learning-based load balancing."""

    def predict_fragment_cost(self, fragment: FragmentTask) -> float:
        """Predict computational cost using trained model."""

    def optimize_distribution(self, fragments: List[FragmentTask]) -> Dict[int, List[FragmentTask]]:
        """Optimize fragment distribution using ML predictions."""
```