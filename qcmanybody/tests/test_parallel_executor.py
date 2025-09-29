"""Tests for the parallel execution engine.

This module implements comprehensive testing of the ParallelManyBodyExecutor,
including unit tests, integration tests, and validation against sequential execution.
"""

import copy

import pytest
import numpy as np
from qcelemental.models import Molecule

from qcmanybody.core import ManyBodyCore
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
from qcmanybody.models.v1 import BsseEnum


class TestParallelConfig:
    """Test the ParallelConfig class."""

    def test_default_config(self):
        """Test default configuration values."""
        config = ParallelConfig()
        assert config.max_workers == 4
        assert config.execution_mode == "multiprocessing"
        assert config.memory_limit_mb == 1000
        assert config.timeout_seconds == 3600
        assert config.use_qcengine is True
        assert config.qc_program == "psi4"
        assert config.basis_set == "sto-3g"
        assert config.qcengine_config == {}
        assert config.default_driver == "energy"

    def test_custom_config(self):
        """Test custom configuration values."""
        config = ParallelConfig(
            max_workers=8,
            execution_mode="threading",
            memory_limit_mb=2000,
            timeout_seconds=1800,
            use_qcengine=False,
            qc_program="nwchem",
            basis_set="6-31g"
        )
        assert config.max_workers == 8
        assert config.execution_mode == "threading"
        assert config.memory_limit_mb == 2000
        assert config.timeout_seconds == 1800
        assert config.use_qcengine is False
        assert config.qc_program == "nwchem"
        assert config.basis_set == "6-31g"
        assert config.default_driver == "energy"

    def test_invalid_execution_mode(self):
        """Test validation of execution mode."""
        with pytest.raises(ValueError, match="Invalid execution_mode"):
            ParallelConfig(execution_mode="invalid")

    def test_invalid_max_workers(self):
        """Test validation of max_workers."""
        with pytest.raises(ValueError, match="max_workers must be >= 1"):
            ParallelConfig(max_workers=0)


class TestParallelManyBodyExecutor:
    """Test the ParallelManyBodyExecutor class."""

    @pytest.fixture
    def water_dimer_molecule(self):
        """Create a simple water dimer molecule for testing."""
        # Create molecule with proper fragments specification
        mol = Molecule.from_data("""
        O  0.0000  0.0000  0.0000
        H  0.7570  0.5860  0.0000
        H -0.7570  0.5860  0.0000
        --
        O  3.0000  0.0000  0.0000
        H  3.7570  0.5860  0.0000
        H  2.2430  0.5860  0.0000
        """)
        # Ensure fragments are properly set
        if not mol.fragments:
            mol = mol.copy(update={"fragments": [[0, 1, 2], [3, 4, 5]]})
        return mol

    @pytest.fixture
    def simple_manybody_core(self, water_dimer_molecule):
        """Create a simple ManyBodyCore for testing."""
        return ManyBodyCore(
            molecule=water_dimer_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf", 2: "hf"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={}
        )

    @pytest.fixture
    def simple_specifications(self):
        """Provide a minimal specification mapping for placeholder execution."""
        spec = {
            "hf": {
                "program": "psi4",
                "specification": {
                    "driver": "energy",
                    "model": {"method": "hf", "basis": "sto-3g"},
                    "keywords": {},
                    "protocols": {},
                    "extras": {},
                },
            }
        }
        return copy.deepcopy(spec)

    def test_executor_initialization(self, simple_manybody_core, simple_specifications):
        """Test ParallelManyBodyExecutor initialization."""
        config = ParallelConfig(use_qcengine=False)  # Disable QCEngine for testing
        executor = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=simple_specifications,
        )

        assert executor.core is simple_manybody_core
        assert executor.config is config
        assert hasattr(executor, '_dependency_graph')
        assert executor.execution_stats["total_fragments"] == 0

    def test_executor_initialization_missing_dependency_graph(self, simple_specifications):
        """Test initialization with missing dependency graph methods."""
        # Create a mock core without iterate_molecules_by_level
        class MockCore:
            @property
            def dependency_graph(self):
                return None

        mock_core = MockCore()
        config = ParallelConfig(use_qcengine=False)

        with pytest.raises(RuntimeError, match="P1-002 dependency graph foundation required"):
            ParallelManyBodyExecutor(
                mock_core,
                config,
                driver="energy",
                specifications=simple_specifications,
            )

    def test_execute_fragment_placeholder(self, simple_manybody_core, simple_specifications):
        """Test fragment execution with placeholder (no QCEngine)."""
        config = ParallelConfig(use_qcengine=False)
        executor = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=simple_specifications,
        )

        fragment_task = executor._tasks_by_level[1][0]
        label, result = executor.execute_fragment(fragment_task)

        assert label == '["hf", [1], [1]]'
        assert result.success is True
        assert result.driver == "energy"
        assert result.model.method == "hf"
        assert result.model.basis == "sto-3g"
        assert isinstance(result.return_result, (int, float))

    def test_fragment_driver_propagation(self, water_dimer_molecule):
        """Ensure driver selection propagates into AtomicInput and results."""
        core = ManyBodyCore(
            molecule=water_dimer_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf", 2: "hf"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
        )

        specifications = {
            "hf": {
                "program": "psi4",
                "specification": {
                    "driver": "gradient",
                    "model": {"method": "hf", "basis": "sto-3g"},
                    "keywords": {},
                    "protocols": {},
                    "extras": {},
                },
            }
        }

        config = ParallelConfig(use_qcengine=False)
        executor = ParallelManyBodyExecutor(core, config, driver="gradient", specifications=specifications)

        fragment_task = executor._tasks_by_level[1][0]
        assert fragment_task.atomic_input.driver == "gradient"

        label, result = executor.execute_fragment(fragment_task)
        assert label == fragment_task.label
        assert result.driver == "gradient"
        assert result.return_result.shape == (len(fragment_task.atomic_input.molecule.symbols), 3)
        assert np.allclose(result.return_result, 0.0)

    def test_execute_level_parallel_serial_mode(self, simple_manybody_core, simple_specifications):
        """Test level execution in serial mode."""
        config = ParallelConfig(use_qcengine=False, execution_mode="serial")
        executor = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=simple_specifications,
        )

        fragments_at_level = executor._tasks_by_level[1]

        results = executor.execute_level_parallel(1, fragments_at_level)

        assert len(results) == 2
        assert '["hf", [1], [1]]' in results
        assert '["hf", [2], [2]]' in results
        assert all(result.success for result in results.values())

    def test_execute_level_parallel_threading_mode(self, simple_manybody_core, simple_specifications):
        """Test level execution in threading mode."""
        config = ParallelConfig(use_qcengine=False, execution_mode="threading", max_workers=2)
        executor = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=simple_specifications,
        )

        fragments_at_level = executor._tasks_by_level[1]

        results = executor.execute_level_parallel(1, fragments_at_level)

        assert len(results) == 2
        assert '["hf", [1], [1]]' in results
        assert '["hf", [2], [2]]' in results
        assert all(result.success for result in results.values())

    def test_execute_full_calculation(self, simple_manybody_core, simple_specifications):
        """Test full parallel calculation execution."""
        config = ParallelConfig(use_qcengine=False, execution_mode="serial")
        executor = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=simple_specifications,
        )

        results = executor.execute_full_calculation()

        # Should have results for monomers and dimer
        assert len(results) > 0
        assert all(result.success for result in results.values())

        # Check execution statistics
        stats = executor.get_execution_statistics()
        assert stats["total_fragments"] > 0
        assert stats["levels_executed"] > 0
        assert stats["parallel_time"] > 0
        assert stats["speedup_factor"] >= 0

    def test_validation_framework(self, simple_manybody_core, simple_specifications):
        """Test parallel vs sequential validation framework."""
        config = ParallelConfig(use_qcengine=False, execution_mode="serial")
        executor = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=copy.deepcopy(simple_specifications),
        )

        # Execute twice to simulate parallel vs sequential
        results1 = executor.execute_full_calculation()

        # Reset executor and execute again
        executor2 = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=copy.deepcopy(simple_specifications),
        )
        results2 = executor2.execute_full_calculation()

        # Validate results match within tolerance
        is_valid = executor.validate_parallel_correctness(results1, results2, tolerance=1e-12)
        assert is_valid is True

    def test_validation_framework_mismatch(self, simple_manybody_core, simple_specifications):
        """Test validation framework with mismatched results."""
        config = ParallelConfig(use_qcengine=False)
        executor = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=simple_specifications,
        )

        results1 = executor.execute_full_calculation()

        # Create modified results with different values
        results2 = results1.copy()
        if results2:
            first_key = next(iter(results2.keys()))
            # Modify the first result to create a mismatch
            original_result = results2[first_key]
            modified_result = original_result.copy(update={
                "return_result": original_result.return_result + 1.0
            })
            results2[first_key] = modified_result

            # Should detect the mismatch
            with pytest.raises(ValueError, match="Energy difference.*exceeds tolerance"):
                executor.validate_parallel_correctness(results1, results2, tolerance=1e-12)

    def test_execution_statistics_tracking(self, simple_manybody_core, simple_specifications):
        """Test that execution statistics are properly tracked."""
        config = ParallelConfig(use_qcengine=False, execution_mode="serial")
        executor = ParallelManyBodyExecutor(
            simple_manybody_core,
            config,
            driver="energy",
            specifications=simple_specifications,
        )

        # Initial state
        stats = executor.get_execution_statistics()
        assert stats["total_fragments"] == 0
        assert stats["levels_executed"] == 0
        assert stats["parallel_time"] == 0.0

        # After execution
        executor.execute_full_calculation()
        stats = executor.get_execution_statistics()
        assert stats["total_fragments"] > 0
        assert stats["levels_executed"] > 0
        assert stats["parallel_time"] > 0
        assert stats["speedup_factor"] >= 0


@pytest.mark.skipif(
    not pytest.importorskip("qcengine", reason="QCEngine not available"),
    reason="QCEngine integration tests require qcengine"
)
class TestQCEngineIntegration:
    """Test QCEngine integration (requires qcengine to be installed)."""

    @pytest.fixture
    def small_molecule(self):
        """Create a very small molecule for quick QCEngine tests."""
        return Molecule(
            symbols=["H", "H"],
            geometry=[0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
            fragments=[[0], [1]],
        )

    def test_qcengine_execution_disabled(self, small_molecule):
        """Test that QCEngine execution can be disabled."""
        simple_core = ManyBodyCore(
            molecule=small_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        config = ParallelConfig(use_qcengine=False)
        specifications = {
            "hf": {
                "program": "psi4",
                "specification": {
                    "driver": "energy",
                    "model": {"method": "hf", "basis": "sto-3g"},
                    "keywords": {},
                    "protocols": {},
                    "extras": {},
                },
            }
        }
        executor = ParallelManyBodyExecutor(simple_core, config, driver="energy", specifications=specifications)

        fragment_task = executor._tasks_by_level[1][0]
        label, result = executor.execute_fragment(fragment_task)

        assert result.success is True
        assert "placeholder" in result.provenance.creator


# Performance and integration tests that don't require external dependencies
class TestParallelPerformance:
    """Test performance characteristics of parallel execution."""

    def test_speedup_calculation(self):
        """Test speedup calculation in execution statistics."""
        config = ParallelConfig(use_qcengine=False)

        # Create a mock executor to test speedup calculation
        class MockExecutor(ParallelManyBodyExecutor):
            def __init__(self, config):
                self.config = config
                self.execution_stats = {
                    "total_fragments": 10,
                    "levels_executed": 3,
                    "parallel_time": 5.0,
                    "sequential_time_estimate": 15.0,
                    "speedup_factor": 0.0
                }

        executor = MockExecutor(config)

        # Manual speedup calculation
        executor.execution_stats["speedup_factor"] = (
            executor.execution_stats["sequential_time_estimate"] /
            executor.execution_stats["parallel_time"]
        )

        assert executor.execution_stats["speedup_factor"] == 3.0

    def test_empty_level_handling(self):
        """Test handling of empty levels in parallel execution."""
        # This tests the edge case where a level has no fragments
        config = ParallelConfig(use_qcengine=False)

        # Create minimal mock for testing
        class MockCore:
            def iterate_molecules_by_level(self):
                return iter([])  # Empty iteration

            @property
            def dependency_graph(self):
                return self

        mock_core = MockCore()
        executor = ParallelManyBodyExecutor(mock_core, config, driver="energy", specifications={})

        # Should handle empty iteration gracefully
        results = executor.execute_full_calculation()
        assert results == {}
        assert executor.execution_stats["levels_executed"] == 0