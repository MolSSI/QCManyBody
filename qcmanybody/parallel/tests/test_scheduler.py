"""
Tests for task scheduling and load balancing.
"""

import pytest
from typing import List

from qcmanybody.parallel.scheduler import (
    TaskScheduler,
    SchedulingStrategy,
    assign_task_priorities,
    estimate_load_balance,
)
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
def create_tasks(mock_molecule, mock_atomic_input):
    """Factory fixture to create test tasks."""
    def _create(n: int = 5) -> List[ParallelTask]:
        """
        Create n test tasks with varying properties.

        Tasks are created with:
        - Varying n-body levels (1-3)
        - Varying estimated costs (1.0-5.0)
        - Varying priorities (0-20)
        """
        tasks = []
        for i in range(n):
            task = ParallelTask(
                task_id=f"task_{i}",
                chemistry="hf/sto-3g",
                label=f"test_{i}",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
                nbody=(i % 3) + 1,  # Vary n-body level: 1, 2, 3, 1, 2, ...
                estimated_cost=float(i % 5) + 1.0,  # Vary cost: 1, 2, 3, 4, 5, 1, ...
                priority=(i * 4) % 21  # Vary priority: 0, 4, 8, 12, 16, 20, 0, ...
            )
            tasks.append(task)
        return tasks
    return _create


# ============================================================================
# Test SchedulingStrategy
# ============================================================================


class TestSchedulingStrategy:
    """Test SchedulingStrategy configuration."""

    def test_default_strategy(self):
        """Test default strategy values."""
        strategy = SchedulingStrategy()
        assert strategy.name == "priority_first"
        assert strategy.enable_load_balancing is True
        assert strategy.chunk_size is None
        assert strategy.reorder_tasks is True

    def test_custom_strategy(self):
        """Test custom strategy configuration."""
        strategy = SchedulingStrategy(
            name="cost_first",
            enable_load_balancing=False,
            chunk_size=10,
            reorder_tasks=False
        )
        assert strategy.name == "cost_first"
        assert strategy.enable_load_balancing is False
        assert strategy.chunk_size == 10
        assert strategy.reorder_tasks is False


# ============================================================================
# Test TaskScheduler Initialization
# ============================================================================


class TestTaskSchedulerInit:
    """Test TaskScheduler initialization."""

    def test_default_init(self):
        """Test default initialization."""
        scheduler = TaskScheduler()
        assert scheduler.strategy.name == "priority_first"
        assert scheduler.n_workers == 1

    def test_custom_init(self):
        """Test initialization with custom strategy."""
        strategy = SchedulingStrategy(name="nbody_first", chunk_size=20)
        scheduler = TaskScheduler(strategy=strategy, n_workers=4)
        assert scheduler.strategy.name == "nbody_first"
        assert scheduler.strategy.chunk_size == 20
        assert scheduler.n_workers == 4


# ============================================================================
# Test Scheduling Strategies
# ============================================================================


class TestSchedulingStrategies:
    """Test different scheduling strategies."""

    def test_priority_first_scheduling(self, create_tasks):
        """Test priority-first scheduling."""
        tasks = create_tasks(10)
        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(name="priority_first"),
            n_workers=4
        )

        scheduled = scheduler.schedule(tasks)

        # Check that tasks are ordered by priority (highest first)
        priorities = [t.priority for t in scheduled]
        assert priorities == sorted(priorities, reverse=True)

    def test_cost_first_scheduling(self, create_tasks):
        """Test cost-first scheduling."""
        tasks = create_tasks(10)
        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(name="cost_first"),
            n_workers=4
        )

        scheduled = scheduler.schedule(tasks)

        # Check that tasks are ordered by cost (highest first)
        costs = [t.estimated_cost for t in scheduled]
        assert costs == sorted(costs, reverse=True)

    def test_nbody_first_scheduling(self, create_tasks):
        """Test n-body-first scheduling."""
        tasks = create_tasks(10)
        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(name="nbody_first"),
            n_workers=4
        )

        scheduled = scheduler.schedule(tasks)

        # Check that tasks are ordered by n-body level (highest first)
        nbody_levels = [t.nbody for t in scheduled]
        assert nbody_levels == sorted(nbody_levels, reverse=True)

    def test_fifo_scheduling(self, create_tasks):
        """Test FIFO (no reordering) scheduling."""
        tasks = create_tasks(10)
        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(name="fifo"),
            n_workers=4
        )

        scheduled = scheduler.schedule(tasks)

        # Check that order is preserved
        original_ids = [t.task_id for t in tasks]
        scheduled_ids = [t.task_id for t in scheduled]
        assert original_ids == scheduled_ids

    def test_dependency_aware_scheduling(self, create_tasks, mock_molecule, mock_atomic_input):
        """Test dependency-aware scheduling."""
        # Create tasks with dependencies
        tasks = []
        for i in range(5):
            depends_on = [f"task_{i-1}"] if i > 0 else []
            task = ParallelTask(
                task_id=f"task_{i}",
                chemistry="hf/sto-3g",
                label=f"test_{i}",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
                depends_on=depends_on,
                priority=i  # Give later tasks higher priority
            )
            tasks.append(task)

        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(name="dependency_aware"),
            n_workers=4
        )

        scheduled = scheduler.schedule(tasks)

        # Check that dependencies are respected
        # Each task should come after its dependencies
        completed_ids = set()
        for task in scheduled:
            for dep_id in task.depends_on:
                assert dep_id in completed_ids, \
                    f"Task {task.task_id} scheduled before dependency {dep_id}"
            completed_ids.add(task.task_id)

    def test_dependency_cycle_fallback(self, mock_molecule, mock_atomic_input):
        """Test fallback to priority scheduling when cycle detected."""
        # Create circular dependency
        tasks = [
            ParallelTask(
                task_id="task_0",
                chemistry="hf/sto-3g",
                label="test_0",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
                depends_on=["task_2"],
                priority=10
            ),
            ParallelTask(
                task_id="task_1",
                chemistry="hf/sto-3g",
                label="test_1",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
                depends_on=["task_0"],
                priority=20
            ),
            ParallelTask(
                task_id="task_2",
                chemistry="hf/sto-3g",
                label="test_2",
                molecule=mock_molecule,
                atomic_input=mock_atomic_input,
                depends_on=["task_1"],
                priority=5
            ),
        ]

        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(name="dependency_aware"),
            n_workers=4
        )

        # Should fallback to priority scheduling
        scheduled = scheduler.schedule(tasks)

        # Check fallback to priority order
        assert len(scheduled) == 3
        priorities = [t.priority for t in scheduled]
        assert priorities == sorted(priorities, reverse=True)


# ============================================================================
# Test Reordering Control
# ============================================================================


class TestReorderingControl:
    """Test task reordering control."""

    def test_reordering_disabled(self, create_tasks):
        """Test that reordering can be disabled."""
        tasks = create_tasks(10)
        original_ids = [t.task_id for t in tasks]

        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(
                name="priority_first",
                reorder_tasks=False
            ),
            n_workers=4
        )

        scheduled = scheduler.schedule(tasks)
        scheduled_ids = [t.task_id for t in scheduled]

        # Order should be preserved
        assert original_ids == scheduled_ids


# ============================================================================
# Test Task Chunking
# ============================================================================


class TestTaskChunking:
    """Test task chunking for batch submission."""

    def test_create_chunks_with_fixed_size(self, create_tasks):
        """Test chunking with fixed chunk size."""
        tasks = create_tasks(20)
        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(chunk_size=5),
            n_workers=4
        )

        chunks = scheduler.create_chunks(tasks)

        assert len(chunks) == 4  # 20 tasks / 5 per chunk
        for chunk in chunks:
            assert len(chunk) == 5

    def test_create_chunks_auto_size(self, create_tasks):
        """Test chunking with automatic size determination."""
        tasks = create_tasks(30)
        scheduler = TaskScheduler(n_workers=4)

        chunks = scheduler.create_chunks(tasks)

        # Should create multiple chunks
        assert len(chunks) > 1

        # All tasks should be included
        total_tasks = sum(len(chunk) for chunk in chunks)
        assert total_tasks == 30

    def test_create_chunks_empty_tasks(self):
        """Test chunking with no tasks."""
        scheduler = TaskScheduler(n_workers=4)
        chunks = scheduler.create_chunks([])
        assert chunks == []

    def test_create_chunks_fewer_tasks_than_workers(self, create_tasks):
        """Test chunking when tasks < workers."""
        tasks = create_tasks(3)
        scheduler = TaskScheduler(n_workers=10)

        chunks = scheduler.create_chunks(tasks)

        # Should create single chunk with all tasks
        assert len(chunks) == 1
        assert len(chunks[0]) == 3


# ============================================================================
# Test Scheduling Statistics
# ============================================================================


class TestSchedulingStats:
    """Test scheduling statistics collection."""

    def test_get_stats(self, create_tasks):
        """Test retrieving scheduling statistics."""
        tasks = create_tasks(10)
        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(name="priority_first"),
            n_workers=4
        )

        scheduler.schedule(tasks)
        stats = scheduler.get_stats()

        assert "strategy" in stats
        assert stats["strategy"] == "priority_first"
        assert stats["n_tasks"] == 10
        assert "reordered" in stats
        assert "priority_range" in stats
        assert "cost_range" in stats
        assert "avg_priority" in stats
        assert "avg_cost" in stats


# ============================================================================
# Test Priority Assignment
# ============================================================================


class TestPriorityAssignment:
    """Test priority assignment functions."""

    def test_assign_nbody_first_priorities(self, create_tasks):
        """Test n-body-first priority assignment."""
        tasks = create_tasks(10)

        # Reset priorities
        for task in tasks:
            task.priority = 0

        assign_task_priorities(tasks, strategy="nbody_first")

        # Higher n-body should have higher priority
        for task in tasks:
            assert task.priority == task.nbody * 10

    def test_assign_cost_first_priorities(self, create_tasks):
        """Test cost-first priority assignment."""
        tasks = create_tasks(10)

        assign_task_priorities(tasks, strategy="cost_first")

        # Higher cost should have higher priority
        costs = [t.estimated_cost for t in tasks]
        priorities = [t.priority for t in tasks]

        # Check correlation (not exact due to normalization)
        max_cost_idx = costs.index(max(costs))
        max_priority_idx = priorities.index(max(priorities))
        assert max_cost_idx == max_priority_idx

    def test_assign_balanced_priorities(self, create_tasks):
        """Test balanced priority assignment."""
        tasks = create_tasks(10)

        assign_task_priorities(tasks, strategy="balanced")

        # All tasks should have assigned priorities
        for task in tasks:
            assert task.priority >= 0
            assert task.priority <= 100

    def test_assign_priorities_empty_tasks(self):
        """Test priority assignment with no tasks."""
        assign_task_priorities([], strategy="nbody_first")
        # Should not raise exception

    def test_assign_priorities_unknown_strategy(self, create_tasks):
        """Test priority assignment with unknown strategy."""
        tasks = create_tasks(5)
        original_priorities = [t.priority for t in tasks]

        assign_task_priorities(tasks, strategy="unknown")

        # Priorities should remain unchanged
        assert [t.priority for t in tasks] == original_priorities


# ============================================================================
# Test Load Balance Estimation
# ============================================================================


class TestLoadBalanceEstimation:
    """Test load balance estimation."""

    def test_estimate_perfect_balance(self, create_tasks):
        """Test estimation with perfectly balanced tasks."""
        # Create tasks with identical costs
        tasks = create_tasks(8)
        for task in tasks:
            task.estimated_cost = 1.0

        result = estimate_load_balance(tasks, n_workers=4)

        assert "ideal_time" in result
        assert "estimated_time" in result
        assert "efficiency" in result
        assert "imbalance" in result

        # Perfect balance: efficiency should be 1.0
        assert result["efficiency"] == pytest.approx(1.0, abs=0.01)
        assert result["imbalance"] == pytest.approx(1.0, abs=0.01)

    def test_estimate_imbalanced_load(self, create_tasks):
        """Test estimation with imbalanced tasks."""
        # Create tasks with very different costs
        tasks = create_tasks(4)
        tasks[0].estimated_cost = 10.0
        tasks[1].estimated_cost = 1.0
        tasks[2].estimated_cost = 1.0
        tasks[3].estimated_cost = 1.0

        result = estimate_load_balance(tasks, n_workers=2)

        # Should show imbalance
        assert result["imbalance"] > 1.0
        assert result["efficiency"] < 1.0

    def test_estimate_empty_tasks(self):
        """Test estimation with no tasks."""
        result = estimate_load_balance([], n_workers=4)

        assert result["ideal_time"] == 0.0
        assert result["estimated_time"] == 0.0
        assert result["efficiency"] == 1.0
        assert result["imbalance"] == 1.0

    def test_estimate_zero_workers(self, create_tasks):
        """Test estimation with zero workers."""
        tasks = create_tasks(5)
        result = estimate_load_balance(tasks, n_workers=0)

        # Should return sensible defaults
        assert result["efficiency"] == 1.0
        assert result["imbalance"] == 1.0


# ============================================================================
# Integration Tests
# ============================================================================


class TestSchedulerIntegration:
    """Integration tests for scheduler workflow."""

    def test_full_scheduling_workflow(self, create_tasks):
        """Test complete scheduling workflow."""
        # Create tasks
        tasks = create_tasks(20)

        # Assign priorities
        assign_task_priorities(tasks, strategy="nbody_first")

        # Create scheduler
        strategy = SchedulingStrategy(
            name="priority_first",
            chunk_size=5
        )
        scheduler = TaskScheduler(strategy=strategy, n_workers=4)

        # Schedule tasks
        scheduled = scheduler.schedule(tasks)
        assert len(scheduled) == 20

        # Get statistics
        stats = scheduler.get_stats()
        assert stats["n_tasks"] == 20

        # Create chunks
        chunks = scheduler.create_chunks(scheduled)
        assert len(chunks) == 4

        # Estimate load balance
        balance = estimate_load_balance(scheduled, n_workers=4)
        assert "efficiency" in balance

    def test_scheduling_preserves_task_properties(self, create_tasks):
        """Test that scheduling preserves all task properties."""
        tasks = create_tasks(5)
        original_props = [
            (t.task_id, t.chemistry, t.nbody, t.estimated_cost)
            for t in tasks
        ]

        scheduler = TaskScheduler(
            strategy=SchedulingStrategy(name="cost_first"),
            n_workers=4
        )

        scheduled = scheduler.schedule(tasks)
        scheduled_props = [
            (t.task_id, t.chemistry, t.nbody, t.estimated_cost)
            for t in scheduled
        ]

        # All properties should be preserved (order may change)
        assert sorted(original_props) == sorted(scheduled_props)
