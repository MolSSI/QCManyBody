# Integration and Deployment Guide

This guide covers integrating the QCManyBody Parallel Execution system into existing workflows and deploying it in production environments.

## 📋 Integration Documentation

### Production Deployment
- **[Production Setup](production-setup.md)** - Production environment configuration
- **[Environment Management](environment-management.md)** - Conda environments and dependencies
- **[HPC Integration](hpc-integration.md)** - High-performance computing deployment
- **[Cloud Deployment](cloud-deployment.md)** - Cloud platform integration

### Workflow Integration
- **[Existing Workflows](workflow-integration.md)** - Integration with current workflows
- **[CI/CD Integration](cicd-integration.md)** - Continuous integration and deployment
- **[Monitoring and Logging](monitoring-logging.md)** - Production monitoring setup
- **[Error Handling](error-handling.md)** - Robust error handling and recovery

## 🚀 Quick Integration Guide

### Step 1: Install Dependencies

#### Basic Installation
```bash
# Core requirements (already available in QCManyBody)
pip install qcelemental numpy pydantic

# For real quantum chemistry
pip install qcengine

# Install quantum chemistry program (choose one)
conda install -c conda-forge psi4        # Psi4 (recommended)
conda install -c conda-forge nwchem      # NWChem
# or other QCEngine-supported programs
```

#### Production Environment
```bash
# Create dedicated conda environment
conda create -n qcmanybody-parallel python=3.9
conda activate qcmanybody-parallel

# Install complete stack
conda install -c conda-forge psi4 qcengine
pip install qcelemental numpy pydantic
```

### Step 2: Basic Integration

#### Replace Existing ManyBodyCore Usage

**Before (sequential only):**
```python
from qcmanybody import ManyBodyCore, BsseEnum

# Traditional sequential execution
core = ManyBodyCore(
    molecule=molecule,
    bsse_type=[BsseEnum.nocp],
    levels={1: "hf", 2: "mp2"},
    return_total_data=False,
    supersystem_ie_only=False,
    embedding_charges={}
)

# Sequential execution (no parallelization)
results = core.compute()  # Uses ManyBodyComputer internally
```

**After (with parallel execution):**
```python
from qcmanybody import ManyBodyCore, BsseEnum
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig

# Same ManyBodyCore setup
core = ManyBodyCore(
    molecule=molecule,
    bsse_type=[BsseEnum.nocp],
    levels={1: "hf", 2: "mp2"},
    return_total_data=False,
    supersystem_ie_only=False,
    embedding_charges={}
)

# Add parallel execution configuration
config = ParallelConfig(
    max_workers=4,
    execution_mode="threading",
    use_qcengine=True,
    qc_program="psi4"
)

# Parallel execution
executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()
```

### Step 3: Workflow Integration

#### Convenience Function Pattern

```python
def run_parallel_manybody(molecule, method_levels, parallel_config=None, **kwargs):
    """Convenience function for parallel many-body calculations."""

    # Set up ManyBodyCore with provided parameters
    core = ManyBodyCore(
        molecule=molecule,
        levels=method_levels,
        **kwargs
    )

    # Use provided config or create default
    if parallel_config is None:
        parallel_config = ParallelConfig(
            max_workers=4,
            execution_mode="threading",
            use_qcengine=True
        )

    # Execute with parallelization
    executor = ParallelManyBodyExecutor(core, parallel_config)
    return executor.execute_full_calculation()

# Usage in existing workflows
results = run_parallel_manybody(
    molecule=water_trimer,
    method_levels={1: "hf", 2: "mp2", 3: "ccsd(t)"},
    bsse_type=[BsseEnum.cp],
    parallel_config=ParallelConfig(max_workers=8)
)
```

## 🏭 Production Deployment

### Production Configuration Template

```python
def create_production_config(environment="hpc"):
    """Create production-optimized parallel configuration."""

    if environment == "hpc":
        # High-performance computing cluster
        return ParallelConfig(
            max_workers=16,
            execution_mode="multiprocessing",
            use_qcengine=True,
            qc_program="psi4",
            memory_limit_mb=4000,
            timeout_seconds=7200,
            qcengine_config={
                "keywords": {"scf_type": "df", "mp2_type": "df"},
                "protocols": {"stdout": False},
                "task_config": {
                    "scratch_directory": "/scratch/$SLURM_JOB_ID",
                    "ncores": 1
                }
            }
        )

    elif environment == "workstation":
        # High-end workstation
        return ParallelConfig(
            max_workers=8,
            execution_mode="threading",
            use_qcengine=True,
            qc_program="psi4",
            memory_limit_mb=2000,
            timeout_seconds=3600,
            qcengine_config={
                "keywords": {"scf_type": "df"},
                "protocols": {"stdout": False}
            }
        )

    elif environment == "cloud":
        # Cloud computing instance
        return ParallelConfig(
            max_workers=4,
            execution_mode="threading",
            use_qcengine=True,
            qc_program="psi4",
            memory_limit_mb=1500,
            timeout_seconds=1800,
            qcengine_config={
                "protocols": {"stdout": False}
            }
        )

    else:
        # Conservative default
        return ParallelConfig(
            max_workers=2,
            execution_mode="threading",
            use_qcengine=True,
            memory_limit_mb=1000,
            timeout_seconds=1800
        )
```

### Environment Detection

```python
def detect_execution_environment():
    """Automatically detect the execution environment."""

    import os
    import shutil

    # Check for HPC environment indicators
    if any(env_var in os.environ for env_var in ["SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID"]):
        return "hpc"

    # Check for cloud environment indicators
    if any(env_var in os.environ for env_var in ["AWS_EXECUTION_ENV", "GOOGLE_CLOUD_PROJECT"]):
        return "cloud"

    # Check CPU count for workstation vs laptop
    import multiprocessing
    cpu_count = multiprocessing.cpu_count()

    if cpu_count >= 16:
        return "workstation"
    elif cpu_count >= 8:
        return "desktop"
    else:
        return "laptop"

# Usage
environment = detect_execution_environment()
config = create_production_config(environment)
```

### Resource Management

```python
def configure_resources(config: ParallelConfig, total_memory_gb=None, total_cores=None):
    """Configure resources based on system capabilities."""

    import psutil
    import multiprocessing

    # Auto-detect system resources if not provided
    if total_memory_gb is None:
        total_memory_gb = psutil.virtual_memory().total / (1024**3)

    if total_cores is None:
        total_cores = multiprocessing.cpu_count()

    # Calculate optimal settings
    # Reserve 25% of memory for system and other processes
    available_memory_gb = total_memory_gb * 0.75

    # Calculate memory per worker
    memory_per_worker_mb = int((available_memory_gb * 1024) / config.max_workers)

    # Ensure workers don't exceed CPU count
    optimal_workers = min(config.max_workers, total_cores)

    # Update configuration
    optimized_config = ParallelConfig(
        max_workers=optimal_workers,
        execution_mode=config.execution_mode,
        use_qcengine=config.use_qcengine,
        qc_program=config.qc_program,
        memory_limit_mb=memory_per_worker_mb,
        timeout_seconds=config.timeout_seconds,
        qcengine_config=config.qcengine_config
    )

    return optimized_config

# Usage
config = create_production_config("workstation")
optimized_config = configure_resources(config)
```

## 🔧 Advanced Integration Patterns

### Workflow Factory Pattern

```python
class ParallelManyBodyWorkflow:
    """Factory for creating parallel many-body workflows."""

    def __init__(self, default_config=None):
        self.default_config = default_config or ParallelConfig()

    def create_executor(self, molecule, method_levels, **kwargs):
        """Create configured executor for given system."""

        # Extract parallel config if provided
        parallel_config = kwargs.pop('parallel_config', self.default_config)

        # Set up ManyBodyCore
        core = ManyBodyCore(
            molecule=molecule,
            levels=method_levels,
            **kwargs
        )

        # Create executor
        return ParallelManyBodyExecutor(core, parallel_config)

    def execute_with_validation(self, molecule, method_levels, **kwargs):
        """Execute with built-in validation."""

        executor = self.create_executor(molecule, method_levels, **kwargs)

        # Run parallel calculation
        parallel_results = executor.execute_full_calculation()

        # Optional validation against sequential
        if kwargs.get('validate_correctness', False):
            sequential_config = ParallelConfig(execution_mode="serial", max_workers=1)
            sequential_executor = self.create_executor(
                molecule, method_levels, parallel_config=sequential_config, **kwargs
            )
            sequential_results = sequential_executor.execute_full_calculation()

            # Validate
            executor.validate_parallel_correctness(
                parallel_results, sequential_results, tolerance=1e-12
            )

        return parallel_results

# Usage
workflow = ParallelManyBodyWorkflow(
    default_config=create_production_config("workstation")
)

results = workflow.execute_with_validation(
    molecule=water_trimer,
    method_levels={1: "hf", 2: "mp2"},
    validate_correctness=True
)
```

### Error Handling and Recovery

```python
def robust_parallel_execution(executor, max_retries=3, fallback_to_serial=True):
    """Execute with error handling and automatic fallback."""

    for attempt in range(max_retries):
        try:
            # Attempt parallel execution
            results = executor.execute_full_calculation()
            return results

        except TimeoutError as e:
            logger.warning(f"Timeout in attempt {attempt + 1}: {e}")
            if attempt < max_retries - 1:
                # Increase timeout for retry
                executor.config.timeout_seconds *= 2
                continue

        except MemoryError as e:
            logger.warning(f"Memory error in attempt {attempt + 1}: {e}")
            if attempt < max_retries - 1:
                # Reduce worker count and memory usage
                executor.config.max_workers = max(1, executor.config.max_workers // 2)
                executor.config.memory_limit_mb = int(executor.config.memory_limit_mb * 0.7)
                continue

        except Exception as e:
            logger.error(f"Unexpected error in attempt {attempt + 1}: {e}")
            if attempt < max_retries - 1:
                continue

    # All retries failed
    if fallback_to_serial:
        logger.info("Falling back to serial execution")
        fallback_config = ParallelConfig(execution_mode="serial", max_workers=1)
        fallback_executor = ParallelManyBodyExecutor(executor.core, fallback_config)
        return fallback_executor.execute_full_calculation()

    else:
        raise RuntimeError("Parallel execution failed after all retries")

# Usage
try:
    results = robust_parallel_execution(executor)
except RuntimeError as e:
    logger.error(f"Complete execution failure: {e}")
    # Handle complete failure (e.g., notify user, log to monitoring)
```

## 📊 Monitoring and Logging

### Production Logging Configuration

```python
def setup_production_logging(log_level="INFO", log_file="qcmanybody_parallel.log"):
    """Configure production logging for parallel execution."""

    import logging
    import logging.handlers

    # Create logger
    logger = logging.getLogger("qcmanybody.parallel")
    logger.setLevel(getattr(logging, log_level.upper()))

    # Create formatters
    detailed_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s'
    )

    simple_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s'
    )

    # File handler for detailed logging
    file_handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=10*1024*1024, backupCount=5
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(detailed_formatter)

    # Console handler for important messages
    console_handler = logging.StreamHandler()
    console_handler.setLevel(getattr(logging, log_level.upper()))
    console_handler.setFormatter(simple_formatter)

    # Add handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

# Usage
logger = setup_production_logging(log_level="INFO")
```

### Performance Monitoring

```python
def monitor_execution_performance(executor):
    """Monitor and log execution performance metrics."""

    # Get execution statistics
    stats = executor.get_execution_statistics()

    # Log performance metrics
    logger.info(f"Execution completed:")
    logger.info(f"  Total fragments: {stats['total_fragments']}")
    logger.info(f"  Execution time: {stats['parallel_time']:.3f}s")
    logger.info(f"  Estimated speedup: {stats['speedup_factor']:.2f}x")

    # Check for performance issues
    if stats['speedup_factor'] < 1.0:
        logger.warning(f"Performance regression detected: {stats['speedup_factor']:.2f}x speedup")

    if 'memory_usage_mb' in stats:
        logger.info(f"  Peak memory usage: {stats['memory_usage_mb']:.1f} MB")

        # Check memory usage
        if stats['memory_usage_mb'] > executor.config.memory_limit_mb * 0.9:
            logger.warning("High memory usage detected")

    return stats

# Usage
executor = ParallelManyBodyExecutor(core, config)
results = executor.execute_full_calculation()
performance_stats = monitor_execution_performance(executor)
```

## 🌐 Cloud Deployment

### AWS Integration

```python
def configure_for_aws(instance_type="c5.4xlarge"):
    """Configure for AWS EC2 deployment."""

    # Instance-specific configurations
    aws_configs = {
        "c5.xlarge": ParallelConfig(max_workers=2, memory_limit_mb=1500),
        "c5.2xlarge": ParallelConfig(max_workers=4, memory_limit_mb=3000),
        "c5.4xlarge": ParallelConfig(max_workers=8, memory_limit_mb=6000),
        "c5.9xlarge": ParallelConfig(max_workers=16, memory_limit_mb=12000),
    }

    base_config = aws_configs.get(instance_type, ParallelConfig())

    # AWS-specific optimizations
    base_config.qcengine_config = {
        "task_config": {
            "scratch_directory": "/tmp/qc_scratch"
        },
        "protocols": {"stdout": False}
    }

    return base_config
```

### Docker Integration

```dockerfile
# Dockerfile for QCManyBody Parallel Execution
FROM continuumio/miniconda3:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment
RUN conda create -n qcmanybody-parallel python=3.9
RUN echo "conda activate qcmanybody-parallel" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Install quantum chemistry software
RUN conda install -c conda-forge psi4 qcengine

# Install QCManyBody with parallel support
COPY . /app/qcmanybody
WORKDIR /app/qcmanybody
RUN pip install -e .

# Set up runtime environment
ENV PYTHONPATH=/app/qcmanybody
ENV CONDA_DEFAULT_ENV=qcmanybody-parallel

# Default command
CMD ["python", "-c", "from qcmanybody.parallel import ParallelManyBodyExecutor; print('QCManyBody Parallel ready')"]
```

## 📋 Best Practices

### Integration Checklist

#### Pre-Integration Testing
- [ ] Validate mathematical correctness with test systems
- [ ] Benchmark performance on representative calculations
- [ ] Test error handling and recovery mechanisms
- [ ] Verify resource usage patterns
- [ ] Test with actual quantum chemistry programs

#### Production Deployment
- [ ] Configure appropriate logging levels
- [ ] Set up monitoring and alerting
- [ ] Implement backup/fallback strategies
- [ ] Document configuration decisions
- [ ] Establish performance baselines

#### Ongoing Maintenance
- [ ] Monitor performance metrics
- [ ] Update configurations based on workload changes
- [ ] Regular validation testing
- [ ] Keep dependencies updated
- [ ] Review and optimize resource usage

### Common Integration Patterns

#### Pattern 1: Drop-in Replacement
```python
# Minimal change to existing code
# Replace ManyBodyComputer with ParallelManyBodyExecutor
```

#### Pattern 2: Configuration-Driven
```python
# Use configuration files to control parallel execution
# Allows easy switching between serial and parallel modes
```

#### Pattern 3: Workflow Integration
```python
# Integrate into existing computational workflows
# Maintain existing APIs while adding parallel capabilities
```

#### Pattern 4: Service-Oriented
```python
# Deploy as a service for multiple users/applications
# Centralized resource management and monitoring
```

---

This integration guide provides comprehensive information for deploying the QCManyBody Parallel Execution system in production environments while maintaining mathematical correctness and optimal performance.