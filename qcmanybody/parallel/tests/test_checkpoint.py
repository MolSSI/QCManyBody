"""
Tests for checkpointing infrastructure.
"""

import pytest
import json
from pathlib import Path
from typing import List

from qcmanybody.parallel.checkpoint import (
    CheckpointManager,
    CheckpointMetadata,
    create_checkpoint_manager,
)
from qcmanybody.parallel.task import ParallelTask, TaskResult
from qcelemental.models import Molecule, AtomicInput


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def temp_checkpoint_file(tmp_path):
    """Create a temporary checkpoint file path."""
    return tmp_path / "checkpoint.json"


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
def create_tasks(mock_molecule, mock_atomic_input):
    """Factory to create test tasks."""
    def _create(n: int = 5) -> List[ParallelTask]:
        tasks = []
        for i in range(n):
            task = ParallelTask(
                task_id=f"task_{i}",
                chemistry="hf/sto-3g",
                label=f"test_{i}",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
            )
            tasks.append(task)
        return tasks
    return _create


@pytest.fixture
def create_results():
    """Factory to create test results."""
    def _create(n: int = 5, all_success: bool = True) -> List[TaskResult]:
        results = []
        for i in range(n):
            # Determine success for this result
            # If all_success is True, all succeed
            # If all_success is False, use alternating pattern (even=success, odd=failure)
            if all_success is True:
                is_success = True
            elif all_success is False:
                is_success = (i % 2 == 0)  # Even indices succeed
            else:
                is_success = all_success  # Use the value directly

            result = TaskResult(
                task_id=f"task_{i}",
                success=is_success,
                error_message=None if is_success else f"Error in task {i}",
                atomic_result={"energy": -5.0 - i} if is_success else None,
                execution_time=1.0 + i * 0.1,
                worker_id=f"worker_{i % 2}"
            )
            results.append(result)
        return results
    return _create


# ============================================================================
# Test CheckpointManager Initialization
# ============================================================================


class TestCheckpointManagerInit:
    """Test CheckpointManager initialization."""

    def test_init_with_path_str(self, tmp_path):
        """Test initialization with string path."""
        checkpoint_file = str(tmp_path / "checkpoint.json")
        manager = CheckpointManager(checkpoint_file)

        assert manager.checkpoint_file == Path(checkpoint_file)
        assert manager.auto_save is True
        assert manager.save_interval == 10

    def test_init_with_path_object(self, temp_checkpoint_file):
        """Test initialization with Path object."""
        manager = CheckpointManager(temp_checkpoint_file)

        assert manager.checkpoint_file == temp_checkpoint_file
        assert manager.auto_save is True

    def test_init_with_custom_settings(self, temp_checkpoint_file):
        """Test initialization with custom settings."""
        manager = CheckpointManager(
            temp_checkpoint_file,
            auto_save=False,
            save_interval=5
        )

        assert manager.auto_save is False
        assert manager.save_interval == 5


# ============================================================================
# Test Checkpoint Existence
# ============================================================================


class TestCheckpointExistence:
    """Test checkpoint file existence checking."""

    def test_exists_when_file_missing(self, temp_checkpoint_file):
        """Test exists() returns False when file doesn't exist."""
        manager = CheckpointManager(temp_checkpoint_file)
        assert not manager.exists()

    def test_exists_when_file_present(self, temp_checkpoint_file):
        """Test exists() returns True when file exists."""
        # Create empty file
        temp_checkpoint_file.touch()

        manager = CheckpointManager(temp_checkpoint_file)
        assert manager.exists()


# ============================================================================
# Test Saving Results
# ============================================================================


class TestSavingResults:
    """Test saving results to checkpoint."""

    def test_save_single_result(self, temp_checkpoint_file, create_results):
        """Test saving a single result."""
        manager = CheckpointManager(temp_checkpoint_file, auto_save=False)
        results = create_results(1)

        manager.save_result(results[0])

        assert manager.has_result("task_0")
        assert len(manager._results) == 1

    def test_save_multiple_results_individually(
        self, temp_checkpoint_file, create_results
    ):
        """Test saving multiple results one at a time."""
        manager = CheckpointManager(temp_checkpoint_file, auto_save=False)
        results = create_results(5)

        for result in results:
            manager.save_result(result)

        assert len(manager._results) == 5
        for i in range(5):
            assert manager.has_result(f"task_{i}")

    def test_save_multiple_results_batch(self, temp_checkpoint_file, create_results):
        """Test saving multiple results in a batch."""
        manager = CheckpointManager(temp_checkpoint_file, auto_save=False)
        results = create_results(5)

        manager.save_results(results)

        assert len(manager._results) == 5

    def test_auto_save_after_interval(self, temp_checkpoint_file, create_results):
        """Test auto-save triggers after save_interval."""
        manager = CheckpointManager(
            temp_checkpoint_file,
            auto_save=True,
            save_interval=3
        )
        results = create_results(5)

        # Save 2 results (below interval)
        manager.save_result(results[0])
        manager.save_result(results[1])
        assert not temp_checkpoint_file.exists()

        # Save 3rd result (hits interval)
        manager.save_result(results[2])
        assert temp_checkpoint_file.exists()

    def test_manual_save(self, temp_checkpoint_file, create_results):
        """Test manual save operation."""
        manager = CheckpointManager(temp_checkpoint_file, auto_save=False)
        results = create_results(3)

        for result in results:
            manager.save_result(result)

        # File shouldn't exist yet
        assert not temp_checkpoint_file.exists()

        # Manual save
        manager.save(total_tasks=3)

        # Now it should exist
        assert temp_checkpoint_file.exists()


# ============================================================================
# Test Loading Results
# ============================================================================


class TestLoadingResults:
    """Test loading results from checkpoint."""

    def test_load_from_nonexistent_file(self, temp_checkpoint_file):
        """Test loading raises error when file doesn't exist."""
        manager = CheckpointManager(temp_checkpoint_file)

        with pytest.raises(FileNotFoundError):
            manager.load()

    def test_load_saved_results(self, temp_checkpoint_file, create_results):
        """Test loading previously saved results."""
        # Save results
        manager1 = CheckpointManager(temp_checkpoint_file)
        results = create_results(5)
        manager1.save_results(results)

        # Load results in new manager
        manager2 = CheckpointManager(temp_checkpoint_file)
        loaded_results = manager2.load()

        assert len(loaded_results) == 5
        for i in range(5):
            assert loaded_results[i].task_id == f"task_{i}"
            assert loaded_results[i].success
            assert loaded_results[i].atomic_result == {"energy": -5.0 - i}

    def test_load_preserves_all_fields(self, temp_checkpoint_file, create_results):
        """Test that all result fields are preserved during save/load."""
        manager1 = CheckpointManager(temp_checkpoint_file)
        original = create_results(1)[0]

        manager1.save_result(original)
        manager1.save()

        manager2 = CheckpointManager(temp_checkpoint_file)
        loaded = manager2.load()[0]

        assert loaded.task_id == original.task_id
        assert loaded.success == original.success
        assert loaded.error_message == original.error_message
        assert loaded.atomic_result == original.atomic_result
        assert loaded.execution_time == original.execution_time
        assert loaded.worker_id == original.worker_id

    def test_load_with_metadata(self, temp_checkpoint_file, create_results):
        """Test loading checkpoint with metadata."""
        manager1 = CheckpointManager(temp_checkpoint_file)
        results = create_results(10)
        manager1.save_results(results)

        manager2 = CheckpointManager(temp_checkpoint_file)
        loaded_results = manager2.load()

        # Check metadata was loaded
        assert manager2._metadata is not None
        assert manager2._metadata.total_tasks == 10
        assert manager2._metadata.completed_tasks == 10
        assert manager2._metadata.failed_tasks == 0


# ============================================================================
# Test Checkpoint File Format
# ============================================================================


class TestCheckpointFileFormat:
    """Test checkpoint file format and structure."""

    def test_checkpoint_file_is_valid_json(self, temp_checkpoint_file, create_results):
        """Test that checkpoint file is valid JSON."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(3)
        manager.save_results(results)

        # Read and parse JSON
        with open(temp_checkpoint_file, 'r') as f:
            data = json.load(f)

        assert "metadata" in data
        assert "results" in data

    def test_checkpoint_metadata_structure(self, temp_checkpoint_file, create_results):
        """Test checkpoint metadata structure."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(5)
        manager.save_results(results)

        with open(temp_checkpoint_file, 'r') as f:
            data = json.load(f)

        metadata = data["metadata"]
        assert "created_at" in metadata
        assert "qcmanybody_version" in metadata
        assert "total_tasks" in metadata
        assert "completed_tasks" in metadata
        assert "failed_tasks" in metadata
        assert "in_progress_tasks" in metadata

        assert metadata["total_tasks"] == 5
        assert metadata["completed_tasks"] == 5

    def test_checkpoint_results_structure(self, temp_checkpoint_file, create_results):
        """Test checkpoint results structure."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(2)
        manager.save_results(results)

        with open(temp_checkpoint_file, 'r') as f:
            data = json.load(f)

        results_dict = data["results"]
        assert len(results_dict) == 2
        assert "task_0" in results_dict
        assert "task_1" in results_dict

        result = results_dict["task_0"]
        assert "task_id" in result
        assert "success" in result
        assert "error_message" in result
        assert "atomic_result" in result


# ============================================================================
# Test Querying Results
# ============================================================================


class TestQueryingResults:
    """Test querying checkpoint for results."""

    def test_get_completed_task_ids(self, temp_checkpoint_file, create_results):
        """Test getting IDs of completed tasks."""
        manager = CheckpointManager(temp_checkpoint_file)
        # Create mixed success/failure results
        results = create_results(5, all_success=False)  # Alternating success/failure
        manager.save_results(results)

        completed_ids = manager.get_completed_task_ids()

        # Should have task_0, task_2, task_4 (even indices)
        assert "task_0" in completed_ids
        assert "task_2" in completed_ids
        assert "task_4" in completed_ids
        assert "task_1" not in completed_ids
        assert "task_3" not in completed_ids

    def test_get_failed_task_ids(self, temp_checkpoint_file, create_results):
        """Test getting IDs of failed tasks."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(5, all_success=False)
        manager.save_results(results)

        failed_ids = manager.get_failed_task_ids()

        # Should have task_1, task_3 (odd indices)
        assert "task_1" in failed_ids
        assert "task_3" in failed_ids
        assert "task_0" not in failed_ids
        assert "task_2" not in failed_ids

    def test_get_result(self, temp_checkpoint_file, create_results):
        """Test getting specific result."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(3)
        manager.save_results(results)

        result = manager.get_result("task_1")

        assert result is not None
        assert result.task_id == "task_1"

    def test_get_nonexistent_result(self, temp_checkpoint_file):
        """Test getting result that doesn't exist."""
        manager = CheckpointManager(temp_checkpoint_file)
        result = manager.get_result("nonexistent")

        assert result is None

    def test_has_result(self, temp_checkpoint_file, create_results):
        """Test checking if result exists."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(2)
        manager.save_results(results)

        assert manager.has_result("task_0")
        assert manager.has_result("task_1")
        assert not manager.has_result("task_2")


# ============================================================================
# Test Task Filtering
# ============================================================================


class TestTaskFiltering:
    """Test filtering pending tasks."""

    def test_filter_pending_tasks(self, temp_checkpoint_file, create_tasks, create_results):
        """Test filtering out completed tasks."""
        manager = CheckpointManager(temp_checkpoint_file)

        # Save results for tasks 0, 1, 2
        results = create_results(3)
        manager.save_results(results)

        # Create all tasks (0-4)
        all_tasks = create_tasks(5)

        # Filter pending
        pending = manager.filter_pending_tasks(all_tasks)

        # Should only have tasks 3 and 4
        assert len(pending) == 2
        assert pending[0].task_id == "task_3"
        assert pending[1].task_id == "task_4"

    def test_filter_all_completed(self, temp_checkpoint_file, create_tasks, create_results):
        """Test filtering when all tasks are completed."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(5)
        manager.save_results(results)

        tasks = create_tasks(5)
        pending = manager.filter_pending_tasks(tasks)

        assert len(pending) == 0

    def test_filter_none_completed(self, temp_checkpoint_file, create_tasks):
        """Test filtering when no tasks are completed."""
        manager = CheckpointManager(temp_checkpoint_file)
        tasks = create_tasks(5)

        pending = manager.filter_pending_tasks(tasks)

        # All tasks should be pending
        assert len(pending) == 5


# ============================================================================
# Test Finalization and Cleanup
# ============================================================================


class TestFinalizationCleanup:
    """Test checkpoint finalization and cleanup."""

    def test_finalize(self, temp_checkpoint_file, create_results):
        """Test finalizing checkpoint."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(3)

        for result in results:
            manager.save_result(result)

        # Finalize (should save and return all results)
        finalized = manager.finalize()

        assert len(finalized) == 3
        assert temp_checkpoint_file.exists()

    def test_clear(self, temp_checkpoint_file, create_results):
        """Test clearing checkpoint data."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(3)
        manager.save_results(results)

        assert len(manager._results) == 3

        manager.clear()

        assert len(manager._results) == 0
        assert manager._save_counter == 0

    def test_delete(self, temp_checkpoint_file, create_results):
        """Test deleting checkpoint file."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(3)
        manager.save_results(results)

        assert temp_checkpoint_file.exists()

        manager.delete()

        assert not temp_checkpoint_file.exists()
        assert len(manager._results) == 0


# ============================================================================
# Test Atomic Writes
# ============================================================================


class TestAtomicWrites:
    """Test atomic write operations."""

    def test_save_uses_temp_file(self, temp_checkpoint_file, create_results):
        """Test that save uses temporary file for atomic write."""
        manager = CheckpointManager(temp_checkpoint_file)
        results = create_results(1)
        manager.save_results(results)

        # Temp file should not exist after successful save
        temp_file = temp_checkpoint_file.with_suffix('.tmp')
        assert not temp_file.exists()
        assert temp_checkpoint_file.exists()


# ============================================================================
# Test Resume Workflow
# ============================================================================


class TestResumeWorkflow:
    """Test checkpoint resume workflow."""

    def test_resume_from_checkpoint(
        self, temp_checkpoint_file, create_tasks, create_results
    ):
        """Test complete resume workflow."""
        # First run: complete some tasks
        manager1 = CheckpointManager(temp_checkpoint_file)
        results_part1 = create_results(3)  # Complete tasks 0, 1, 2
        manager1.save_results(results_part1)

        # Second run: resume
        manager2 = CheckpointManager(temp_checkpoint_file)
        manager2.load()

        all_tasks = create_tasks(5)
        pending_tasks = manager2.filter_pending_tasks(all_tasks)

        # Should only need to run tasks 3 and 4
        assert len(pending_tasks) == 2
        assert pending_tasks[0].task_id == "task_3"
        assert pending_tasks[1].task_id == "task_4"

        # Complete remaining tasks
        results_part2 = create_results(2)
        results_part2[0].task_id = "task_3"
        results_part2[1].task_id = "task_4"
        manager2.save_results(results_part2)

        # All tasks should be complete now
        final_results = manager2.finalize()
        assert len(final_results) == 5


# ============================================================================
# Test Convenience Functions
# ============================================================================


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_create_checkpoint_manager_no_resume(self, temp_checkpoint_file):
        """Test creating manager without resume."""
        manager = create_checkpoint_manager(
            temp_checkpoint_file,
            resume_if_exists=False
        )

        assert isinstance(manager, CheckpointManager)
        assert manager.checkpoint_file == temp_checkpoint_file

    def test_create_checkpoint_manager_with_resume(
        self, temp_checkpoint_file, create_results
    ):
        """Test creating manager with auto-resume."""
        # Create checkpoint file first
        manager1 = CheckpointManager(temp_checkpoint_file)
        results = create_results(3)
        manager1.save_results(results)

        # Create with resume
        manager2 = create_checkpoint_manager(
            temp_checkpoint_file,
            resume_if_exists=True
        )

        # Should have loaded results
        assert len(manager2._results) == 3

    def test_create_checkpoint_manager_resume_missing_file(self, temp_checkpoint_file):
        """Test creating manager with resume when file doesn't exist."""
        manager = create_checkpoint_manager(
            temp_checkpoint_file,
            resume_if_exists=True
        )

        # Should not raise error, just start fresh
        assert isinstance(manager, CheckpointManager)
        assert len(manager._results) == 0


# ============================================================================
# Error Handling Tests
# ============================================================================


class TestErrorHandling:
    """Test error handling in checkpoint operations."""

    def test_load_corrupted_json(self, temp_checkpoint_file):
        """Test loading checkpoint with corrupted JSON."""
        # Write invalid JSON
        with open(temp_checkpoint_file, 'w') as f:
            f.write("{ invalid json ")

        manager = CheckpointManager(temp_checkpoint_file)

        with pytest.raises(ValueError, match="Failed to load checkpoint"):
            manager.load()

    def test_save_to_nonexistent_directory(self, tmp_path, create_results):
        """Test saving to directory that doesn't exist yet."""
        # Path with nested nonexistent directories
        checkpoint_file = tmp_path / "subdir" / "nested" / "checkpoint.json"

        manager = CheckpointManager(checkpoint_file)
        results = create_results(1)
        manager.save_results(results)

        # Should create directories and save successfully
        assert checkpoint_file.exists()
