"""
Tests for ParallelTask and TaskResult data models.
"""

import pytest
import time
from qcmanybody.parallel import ParallelTask, TaskResult
from qcmanybody.parallel.task import TaskStatus


class TestParallelTask:
    """Tests for ParallelTask dataclass."""

    def test_create_basic_task(self, mock_atomic_input):
        """Test creating a basic ParallelTask."""
        task = ParallelTask(
            task_id="test_1",
            chemistry="hf/sto-3g",
            label='["hf/sto-3g", [1], [1]]',
            molecule=mock_atomic_input.molecule,
            atomic_input=mock_atomic_input,
        )

        assert task.task_id == "test_1"
        assert task.chemistry == "hf/sto-3g"
        assert task.priority == 0
        assert task.estimated_cost == 1.0
        assert task.nbody == 1
        assert task.depends_on == []

    def test_create_task_with_metadata(self, mock_atomic_input):
        """Test creating task with all metadata."""
        task = ParallelTask(
            task_id="test_2",
            chemistry="mp2/cc-pvdz",
            label='["mp2/cc-pvdz", [1, 2], [1, 2]]',
            molecule=mock_atomic_input.molecule,
            atomic_input=mock_atomic_input,
            priority=10,
            estimated_cost=5.0,
            nbody=2,
            depends_on=["test_1"],
            metadata={"extra": "info"},
        )

        assert task.priority == 10
        assert task.estimated_cost == 5.0
        assert task.nbody == 2
        assert task.depends_on == ["test_1"]
        assert task.metadata == {"extra": "info"}

    def test_task_hash(self, mock_parallel_task):
        """Test that tasks can be hashed by ID."""
        task1 = mock_parallel_task
        task2 = ParallelTask(
            task_id=task1.task_id,  # Same ID
            chemistry=task1.chemistry,
            label=task1.label,
            molecule=task1.molecule,
            atomic_input=task1.atomic_input,
        )

        assert hash(task1) == hash(task2)

        # Can be used in sets
        task_set = {task1, task2}
        assert len(task_set) == 1

    def test_task_equality(self, mock_parallel_task):
        """Test task equality by ID."""
        task1 = mock_parallel_task
        task2 = ParallelTask(
            task_id=task1.task_id,  # Same ID
            chemistry=task1.chemistry,
            label=task1.label,
            molecule=task1.molecule,
            atomic_input=task1.atomic_input,
        )

        assert task1 == task2

    def test_task_comparison_for_priority(self, mock_atomic_input):
        """Test task comparison for priority queue."""
        task_low = ParallelTask(
            task_id="low",
            chemistry="hf/sto-3g",
            label="low",
            molecule=mock_atomic_input.molecule,
            atomic_input=mock_atomic_input,
            priority=1,
        )

        task_high = ParallelTask(
            task_id="high",
            chemistry="hf/sto-3g",
            label="high",
            molecule=mock_atomic_input.molecule,
            atomic_input=mock_atomic_input,
            priority=10,
        )

        # Higher priority should be "less than" (comes first in queue)
        assert task_high < task_low

    def test_task_repr(self, mock_parallel_task):
        """Test task string representation."""
        repr_str = repr(mock_parallel_task)

        assert "ParallelTask" in repr_str
        assert mock_parallel_task.task_id in repr_str
        assert mock_parallel_task.chemistry in repr_str


class TestTaskResult:
    """Tests for TaskResult dataclass."""

    def test_create_successful_result(self):
        """Test creating a successful TaskResult."""
        result = TaskResult(
            task_id="test_1",
            success=True,
            status=TaskStatus.COMPLETED,
            execution_time=10.5,
            worker_id="worker_0",
        )

        assert result.task_id == "test_1"
        assert result.success
        assert result.status == TaskStatus.COMPLETED
        assert result.execution_time == 10.5
        assert result.worker_id == "worker_0"
        assert result.error_type is None
        assert result.error_message is None

    def test_create_failed_result(self):
        """Test creating a failed TaskResult."""
        result = TaskResult(
            task_id="test_2",
            success=False,
            status=TaskStatus.FAILED,
            error_type="ValueError",
            error_message="Invalid input",
            execution_time=1.0,
        )

        assert result.task_id == "test_2"
        assert not result.success
        assert result.status == TaskStatus.FAILED
        assert result.error_type == "ValueError"
        assert result.error_message == "Invalid input"
        assert result.atomic_result is None

    def test_create_timeout_result(self):
        """Test creating a timeout TaskResult."""
        result = TaskResult(
            task_id="test_3",
            success=False,
            status=TaskStatus.TIMEOUT,
            error_type="TimeoutError",
            error_message="Task exceeded 3600s timeout",
        )

        assert result.status == TaskStatus.TIMEOUT
        assert not result.success

    def test_return_result_property(self):
        """Test return_result property extraction."""
        # Mock atomic result with return_result attribute
        class MockAtomicResult:
            def __init__(self):
                self.return_result = -5.5

        result = TaskResult(
            task_id="test_4",
            success=True,
            atomic_result=MockAtomicResult(),
        )

        assert result.return_result == -5.5

    def test_return_result_none_when_failed(self):
        """Test return_result is None for failed tasks."""
        result = TaskResult(task_id="test_5", success=False)

        assert result.return_result is None

    def test_total_time_property(self):
        """Test total_time property."""
        result = TaskResult(
            task_id="test_6",
            success=True,
            execution_time=10.0,
            queue_time=5.0,
        )

        assert result.total_time == 15.0

    def test_to_dict(self):
        """Test conversion to dictionary."""
        result = TaskResult(
            task_id="test_7",
            success=True,
            status=TaskStatus.COMPLETED,
            execution_time=10.5,
            queue_time=2.5,
            worker_id="worker_1",
            attempt_number=1,
        )

        result_dict = result.to_dict()

        assert isinstance(result_dict, dict)
        assert result_dict["task_id"] == "test_7"
        assert result_dict["success"] is True
        assert result_dict["status"] == "completed"
        assert result_dict["execution_time"] == 10.5
        assert result_dict["queue_time"] == 2.5
        assert result_dict["worker_id"] == "worker_1"

    def test_result_repr(self):
        """Test result string representation."""
        result = TaskResult(
            task_id="test_8",
            success=True,
            status=TaskStatus.COMPLETED,
            execution_time=10.5,
        )

        repr_str = repr(result)

        assert "TaskResult" in repr_str
        assert "test_8" in repr_str
        assert "completed" in repr_str
