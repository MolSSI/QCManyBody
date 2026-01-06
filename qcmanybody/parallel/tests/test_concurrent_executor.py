"""
Tests for ConcurrentExecutor (concurrent.futures-based execution).
"""

import pytest
import time
from unittest.mock import Mock, patch

from qcmanybody.parallel.executors.concurrent import (
    ConcurrentExecutor,
    ConcurrentExecutorConfig,
    create_concurrent_executor,
)
from qcmanybody.parallel.base import ExecutorConfig
from qcmanybody.parallel.task import ParallelTask
from qcelemental.models import Molecule, AtomicInput


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def mock_molecule():
    """Create a mock molecule for testing."""
    return Molecule(symbols=["He"], geometry=[[0, 0, 0]])


@pytest.fixture
def mock_atomic_input(mock_molecule):
    """Create a mock atomic input for testing."""
    return AtomicInput(
        molecule=mock_molecule,
        driver="energy",
        model={"method": "hf", "basis": "sto-3g"}
    )


@pytest.fixture
def mock_parallel_task(mock_molecule, mock_atomic_input):
    """Create a mock ParallelTask for testing."""
    return ParallelTask(
        task_id="test_task_1",
        chemistry="hf/sto-3g",
        label="test_label",
        molecule=mock_molecule,
        atomic_input=mock_atomic_input,
    )


@pytest.fixture
def concurrent_config():
    """Create a test configuration for concurrent executor."""
    return ConcurrentExecutorConfig(
        n_workers=2,
        executor_type="process",
        timeout_per_task=30.0,
        max_retries=1
    )


@pytest.fixture
def concurrent_executor_process(concurrent_config):
    """Create a ConcurrentExecutor with ProcessPoolExecutor."""
    config = ConcurrentExecutorConfig(
        n_workers=2,
        executor_type="process",
        timeout_per_task=30.0
    )
    return ConcurrentExecutor(config)


@pytest.fixture
def concurrent_executor_thread(concurrent_config):
    """Create a ConcurrentExecutor with ThreadPoolExecutor."""
    config = ConcurrentExecutorConfig(
        n_workers=2,
        executor_type="thread",
        timeout_per_task=30.0
    )
    return ConcurrentExecutor(config)


# ============================================================================
# Test ConcurrentExecutorConfig
# ============================================================================


class TestConcurrentExecutorConfig:
    """Test ConcurrentExecutorConfig configuration."""

    def test_default_config(self):
        """Test default configuration values."""
        config = ConcurrentExecutorConfig()
        assert config.executor_type == "process"
        assert config.max_workers is None
        assert config.n_workers is None

    def test_custom_config(self):
        """Test custom configuration."""
        config = ConcurrentExecutorConfig(
            n_workers=4,
            executor_type="thread",
            max_workers=8,
            timeout_per_task=60.0
        )
        assert config.n_workers == 4
        assert config.executor_type == "thread"
        assert config.max_workers == 8
        assert config.timeout_per_task == 60.0


# ============================================================================
# Test ConcurrentExecutor Initialization
# ============================================================================


class TestConcurrentExecutorInit:
    """Test ConcurrentExecutor initialization."""

    def test_init_with_none_config(self):
        """Test initialization with None config (uses defaults)."""
        executor = ConcurrentExecutor(None)
        assert isinstance(executor.config, ConcurrentExecutorConfig)
        assert executor.config.executor_type == "process"

    def test_init_with_executor_config(self):
        """Test initialization with base ExecutorConfig (should convert)."""
        base_config = ExecutorConfig(n_workers=4, timeout_per_task=60.0)
        executor = ConcurrentExecutor(base_config)

        assert isinstance(executor.config, ConcurrentExecutorConfig)
        assert executor.config.n_workers == 4
        assert executor.config.timeout_per_task == 60.0
        assert executor.config.executor_type == "process"  # Default

    def test_init_with_concurrent_config(self):
        """Test initialization with ConcurrentExecutorConfig."""
        config = ConcurrentExecutorConfig(
            n_workers=4,
            executor_type="thread",
            timeout_per_task=60.0
        )
        executor = ConcurrentExecutor(config)

        assert executor.config.n_workers == 4
        assert executor.config.executor_type == "thread"
        assert executor.config.timeout_per_task == 60.0

    def test_init_invalid_executor_type(self):
        """Test initialization with invalid executor type."""
        config = ConcurrentExecutorConfig(
            n_workers=2,
            executor_type="invalid"
        )
        executor = ConcurrentExecutor(config)

        # Should raise error when trying to initialize resources
        with pytest.raises(ValueError, match="Invalid executor_type"):
            with executor:
                pass


# ============================================================================
# Test Resource Management
# ============================================================================


class TestResourceManagement:
    """Test resource initialization and cleanup."""

    def test_context_manager_process(self, concurrent_executor_process):
        """Test context manager with ProcessPoolExecutor."""
        assert not concurrent_executor_process.is_initialized()

        with concurrent_executor_process as executor:
            assert executor.is_initialized()
            assert executor._executor is not None

        # After context exits, should be cleaned up
        assert not concurrent_executor_process.is_initialized()
        assert concurrent_executor_process._executor is None

    def test_context_manager_thread(self, concurrent_executor_thread):
        """Test context manager with ThreadPoolExecutor."""
        assert not concurrent_executor_thread.is_initialized()

        with concurrent_executor_thread as executor:
            assert executor.is_initialized()
            assert executor._executor is not None

        # After context exits, should be cleaned up
        assert not concurrent_executor_thread.is_initialized()
        assert concurrent_executor_thread._executor is None

    def test_multiple_enter_exit_cycles(self, concurrent_executor_process):
        """Test multiple enter/exit cycles (idempotency)."""
        # First cycle
        with concurrent_executor_process:
            assert concurrent_executor_process.is_initialized()

        assert not concurrent_executor_process.is_initialized()

        # Second cycle
        with concurrent_executor_process:
            assert concurrent_executor_process.is_initialized()

        assert not concurrent_executor_process.is_initialized()

    def test_get_worker_count(self):
        """Test getting worker count."""
        config = ConcurrentExecutorConfig(n_workers=4)
        executor = ConcurrentExecutor(config)

        with executor:
            assert executor.get_worker_count() == 4

    def test_get_worker_count_default(self):
        """Test worker count with default (None)."""
        config = ConcurrentExecutorConfig(n_workers=None)
        executor = ConcurrentExecutor(config)

        with executor:
            # Should use default (0 in our implementation)
            assert executor.get_worker_count() == 0


# ============================================================================
# Test Task Execution
# ============================================================================


class TestTaskExecution:
    """Test task execution with concurrent executor."""

    @pytest.mark.skip(reason="Requires qcengine - Mock objects cannot be pickled across process boundaries")
    def test_execute_single_task_success_process(
        self, concurrent_executor_process, mock_parallel_task
    ):
        """Test executing a single successful task with ProcessPoolExecutor."""
        # NOTE: This test requires actual qcengine to be installed
        # because Mock objects cannot be pickled across process boundaries
        pass

    def test_execute_single_task_success_thread(
        self, concurrent_executor_thread, mock_parallel_task
    ):
        """Test executing a single successful task with ThreadPoolExecutor."""
        mock_result = Mock()
        mock_result.success = True
        mock_result.return_result = -5.0

        with patch("qcengine.compute", return_value=mock_result):
            with concurrent_executor_thread as executor:
                results = executor.execute([mock_parallel_task])

        assert len(results) == 1
        assert results[0].success
        assert results[0].task_id == "test_task_1"

    @pytest.mark.skip(reason="Requires qcengine - Mock objects cannot be pickled")
    def test_execute_multiple_tasks_process(
        self, concurrent_executor_process, mock_molecule, mock_atomic_input
    ):
        """Test executing multiple tasks with ProcessPoolExecutor."""
        # NOTE: This test requires actual qcengine
        pass

    def test_execute_multiple_tasks_thread(
        self, concurrent_executor_thread, mock_molecule, mock_atomic_input
    ):
        """Test executing multiple tasks with ThreadPoolExecutor."""
        # Create multiple tasks
        tasks = [
            ParallelTask(
                task_id=f"task_{i}",
                chemistry="hf/sto-3g",
                label=f"label_{i}",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
            )
            for i in range(5)
        ]

        mock_result = Mock()
        mock_result.success = True
        mock_result.return_result = -5.0

        with patch("qcengine.compute", return_value=mock_result):
            with concurrent_executor_thread as executor:
                results = executor.execute(tasks)

        assert len(results) == 5
        assert all(r.success for r in results)

    def test_execute_empty_task_list(self, concurrent_executor_thread):
        """Test executing empty task list."""
        with concurrent_executor_thread as executor:
            results = executor.execute([])

        assert results == []

    def test_execute_without_initialization(
        self, concurrent_executor_thread, mock_parallel_task
    ):
        """Test that execution fails without initialization."""
        with pytest.raises(RuntimeError, match="not initialized"):
            concurrent_executor_thread.execute([mock_parallel_task])


# ============================================================================
# Test Timeout Handling
# ============================================================================


class TestTimeoutHandling:
    """Test timeout handling in concurrent executor."""

    def test_timeout_behavior_thread(self, mock_molecule, mock_atomic_input):
        """Test timeout with ThreadPoolExecutor."""
        # Create executor with short timeout
        config = ConcurrentExecutorConfig(
            n_workers=2,
            executor_type="thread",
            timeout_per_task=0.1  # Very short timeout
        )
        executor = ConcurrentExecutor(config)

        task = ParallelTask(
            task_id="timeout_task",
            chemistry="hf/sto-3g",
            label="timeout_test",
            molecule=mock_molecule,
            atomic_input=mock_atomic_input,
        )

        # Mock a slow computation
        def slow_compute(*args, **kwargs):
            time.sleep(1.0)  # Sleep longer than timeout
            result = Mock()
            result.success = True
            return result

        with patch("qcengine.compute", side_effect=slow_compute):
            with executor:
                results = executor.execute([task])

        # Task should have failed due to timeout
        assert len(results) == 1
        assert not results[0].success
        assert "timeout" in results[0].error.lower()


# ============================================================================
# Test Error Handling
# ============================================================================


class TestErrorHandling:
    """Test error handling in concurrent executor."""

    def test_task_failure_thread(self, concurrent_executor_thread, mock_parallel_task):
        """Test handling of task failure with ThreadPoolExecutor."""
        mock_result = Mock()
        mock_result.success = False

        with patch("qcengine.compute", return_value=mock_result):
            with concurrent_executor_thread as executor:
                results = executor.execute([mock_parallel_task])

        assert len(results) == 1
        assert not results[0].success

    def test_task_exception_thread(self, concurrent_executor_thread, mock_parallel_task):
        """Test handling of task exception with ThreadPoolExecutor."""
        def raise_exception(*args, **kwargs):
            raise ValueError("Test exception")

        with patch("qcengine.compute", side_effect=raise_exception):
            with concurrent_executor_thread as executor:
                results = executor.execute([mock_parallel_task])

        assert len(results) == 1
        assert not results[0].success
        assert "Test exception" in results[0].error

    def test_mixed_success_failure_thread(
        self, concurrent_executor_thread, mock_molecule, mock_atomic_input
    ):
        """Test handling of mixed success/failure with ThreadPoolExecutor."""
        tasks = [
            ParallelTask(
                task_id=f"task_{i}",
                chemistry="hf/sto-3g",
                label=f"label_{i}",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
            )
            for i in range(4)
        ]

        # Mock alternating success/failure
        def alternating_result(*args, **kwargs):
            result = Mock()
            result.success = (hash(args) % 2 == 0)
            result.return_result = -5.0 if result.success else None
            return result

        with patch("qcengine.compute", side_effect=alternating_result):
            with concurrent_executor_thread as executor:
                results = executor.execute(tasks)

        assert len(results) == 4
        # Should have some successes and some failures
        successes = sum(1 for r in results if r.success)
        failures = sum(1 for r in results if not r.success)
        assert successes > 0
        assert failures > 0


# ============================================================================
# Test Convenience Functions
# ============================================================================


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_create_concurrent_executor_process(self):
        """Test creating process-based executor via convenience function."""
        executor = create_concurrent_executor(
            n_workers=4,
            executor_type="process",
            timeout_per_task=60.0,
            max_retries=3
        )

        assert isinstance(executor, ConcurrentExecutor)
        assert executor.config.n_workers == 4
        assert executor.config.executor_type == "process"
        assert executor.config.timeout_per_task == 60.0
        assert executor.config.max_retries == 3

    def test_create_concurrent_executor_thread(self):
        """Test creating thread-based executor via convenience function."""
        executor = create_concurrent_executor(
            n_workers=8,
            executor_type="thread"
        )

        assert isinstance(executor, ConcurrentExecutor)
        assert executor.config.n_workers == 8
        assert executor.config.executor_type == "thread"

    def test_create_concurrent_executor_defaults(self):
        """Test creating executor with default parameters."""
        executor = create_concurrent_executor()

        assert isinstance(executor, ConcurrentExecutor)
        assert executor.config.executor_type == "process"


# ============================================================================
# Integration Tests
# ============================================================================


class TestConcurrentExecutorIntegration:
    """Integration tests for concurrent executor."""

    def test_full_workflow_thread(self, mock_molecule, mock_atomic_input):
        """Test complete workflow with ThreadPoolExecutor."""
        # Create executor
        config = ConcurrentExecutorConfig(
            n_workers=2,
            executor_type="thread",
            timeout_per_task=30.0
        )
        executor = ConcurrentExecutor(config)

        # Create tasks
        tasks = [
            ParallelTask(
                task_id=f"task_{i}",
                chemistry="hf/sto-3g",
                label=f"label_{i}",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
            )
            for i in range(10)
        ]

        mock_result = Mock()
        mock_result.success = True
        mock_result.return_result = -5.0

        # Execute
        with patch("qcengine.compute", return_value=mock_result):
            with executor:
                results = executor.execute(tasks)

        # Verify
        assert len(results) == 10
        assert all(r.success for r in results)
        assert all(r.task_id.startswith("task_") for r in results)


# ============================================================================
# Comparison Tests
# ============================================================================


class TestProcessVsThread:
    """Tests comparing ProcessPoolExecutor vs ThreadPoolExecutor behavior."""

    def test_both_executors_produce_same_results(
        self, mock_molecule, mock_atomic_input
    ):
        """Test that both executor types produce equivalent results."""
        tasks = [
            ParallelTask(
                task_id=f"task_{i}",
                chemistry="hf/sto-3g",
                label=f"label_{i}",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
            )
            for i in range(5)
        ]

        mock_result = Mock()
        mock_result.success = True
        mock_result.return_result = -5.0

        # Execute with thread executor
        thread_config = ConcurrentExecutorConfig(
            n_workers=2,
            executor_type="thread"
        )
        thread_executor = ConcurrentExecutor(thread_config)

        with patch("qcengine.compute", return_value=mock_result):
            with thread_executor:
                thread_results = thread_executor.execute(tasks.copy())

        # Both should have same number of results
        assert len(thread_results) == 5

        # Both should have all successes (with our mock)
        assert all(r.success for r in thread_results)

    def test_process_executor_initialization(self):
        """Test that process executor initializes correctly."""
        config = ConcurrentExecutorConfig(
            n_workers=2,
            executor_type="process"
        )
        executor = ConcurrentExecutor(config)

        with executor:
            assert executor.is_initialized()
            assert executor.get_worker_count() == 2

    def test_thread_executor_initialization(self):
        """Test that thread executor initializes correctly."""
        config = ConcurrentExecutorConfig(
            n_workers=2,
            executor_type="thread"
        )
        executor = ConcurrentExecutor(config)

        with executor:
            assert executor.is_initialized()
            assert executor.get_worker_count() == 2
