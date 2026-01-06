"""
Tests for base parallel execution classes and interfaces.
"""

import pytest
from qcmanybody.parallel import BaseParallelExecutor, ExecutorConfig


class TestExecutorConfig:
    """Tests for ExecutorConfig dataclass."""

    def test_default_config(self):
        """Test default configuration values."""
        config = ExecutorConfig()

        assert config.n_workers is None
        assert config.timeout_per_task == 3600.0
        assert config.max_retries == 2
        assert config.checkpoint_interval == 10
        assert config.checkpoint_file is None
        assert config.cache_dir is None
        assert config.log_level == "INFO"

    def test_custom_config(self):
        """Test custom configuration values."""
        config = ExecutorConfig(
            n_workers=4,
            timeout_per_task=1800.0,
            max_retries=3,
            checkpoint_interval=20,
            log_level="DEBUG",
        )

        assert config.n_workers == 4
        assert config.timeout_per_task == 1800.0
        assert config.max_retries == 3
        assert config.checkpoint_interval == 20
        assert config.log_level == "DEBUG"

    def test_invalid_n_workers(self):
        """Test that invalid n_workers raises error."""
        with pytest.raises(ValueError, match="n_workers must be >= 1"):
            ExecutorConfig(n_workers=0)

        with pytest.raises(ValueError, match="n_workers must be >= 1"):
            ExecutorConfig(n_workers=-1)

    def test_invalid_timeout(self):
        """Test that invalid timeout raises error."""
        with pytest.raises(ValueError, match="timeout_per_task must be > 0"):
            ExecutorConfig(timeout_per_task=0)

        with pytest.raises(ValueError, match="timeout_per_task must be > 0"):
            ExecutorConfig(timeout_per_task=-1)

    def test_invalid_max_retries(self):
        """Test that invalid max_retries raises error."""
        with pytest.raises(ValueError, match="max_retries must be >= 0"):
            ExecutorConfig(max_retries=-1)

    def test_invalid_checkpoint_interval(self):
        """Test that invalid checkpoint_interval raises error."""
        with pytest.raises(ValueError, match="checkpoint_interval must be >= 1"):
            ExecutorConfig(checkpoint_interval=0)


class TestBaseParallelExecutor:
    """Tests for BaseParallelExecutor interface."""

    def test_cannot_instantiate_abstract_class(self):
        """Test that BaseParallelExecutor cannot be instantiated directly."""
        with pytest.raises(TypeError):
            BaseParallelExecutor()

    def test_executor_name(self, sequential_executor):
        """Test executor name property."""
        assert sequential_executor.name == "SequentialExecutor"

    def test_is_initialized_property(self, sequential_executor):
        """Test is_initialized property."""
        assert not sequential_executor.is_initialized

        sequential_executor.initialize()
        assert sequential_executor.is_initialized

        sequential_executor.shutdown()
        assert not sequential_executor.is_initialized

    def test_context_manager(self, sequential_executor):
        """Test executor as context manager."""
        assert not sequential_executor.is_initialized

        with sequential_executor as ex:
            assert ex.is_initialized
            assert ex is sequential_executor

        assert not sequential_executor.is_initialized

    def test_get_info(self, sequential_executor):
        """Test get_info method."""
        info = sequential_executor.get_info()

        assert isinstance(info, dict)
        assert "name" in info
        assert "n_workers" in info
        assert "is_initialized" in info
        assert "config" in info

        assert info["name"] == "SequentialExecutor"
        assert not info["is_initialized"]

    def test_validate_tasks_empty(self, sequential_executor):
        """Test that empty task list raises error."""
        with pytest.raises(ValueError, match="Task list cannot be empty"):
            sequential_executor.validate_tasks([])

    def test_validate_tasks_duplicate_ids(self, sequential_executor, mock_parallel_task):
        """Test that duplicate task IDs raise error."""
        tasks = [mock_parallel_task, mock_parallel_task]  # Same task twice

        with pytest.raises(ValueError, match="Duplicate task IDs found"):
            sequential_executor.validate_tasks(tasks)

    def test_validate_tasks_valid(self, sequential_executor, mock_parallel_task):
        """Test that valid task list passes validation."""
        from qcmanybody.parallel import ParallelTask

        task2 = ParallelTask(
            task_id="test_task_2",
            chemistry=mock_parallel_task.chemistry,
            label='["hf/sto-3g", [2], [2]]',
            molecule=mock_parallel_task.molecule,
            atomic_input=mock_parallel_task.atomic_input,
        )

        tasks = [mock_parallel_task, task2]

        # Should not raise
        sequential_executor.validate_tasks(tasks)
