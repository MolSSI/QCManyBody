"""
Tests for CLI parallel execution integration.

Tests the CLI commands and argument parsing for parallel execution features.
"""

import pytest
import json
import tempfile
from copy import deepcopy
from pathlib import Path
from argparse import Namespace
from unittest.mock import Mock, patch, MagicMock


# Sample test input data
BASIC_INPUT_PARALLEL = {
    "schema_name": "qcmanybody_cli_input",
    "schema_version": 1,
    "molecule": {
        "source": "inline",
        "inline": {
            "symbols": ["He", "He"],
            "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0]],
            "fragments": [[0], [1]],
            "molecular_charge": 0.0,
            "molecular_multiplicity": 1,
            "units": "angstrom"
        }
    },
    "calculation": {
        "type": "single",
        "single": {
            "driver": "energy",
            "method": "hf",
            "basis": "sto-3g",
            "program": "psi4"
        }
    },
    "bsse": {
        "type": ["cp"]
    },
    "manybody": {
        "max_nbody": 2
    },
    "execution": {
        "parallel": True,
        "n_workers": 2,
        "executor_type": "multiprocessing",
        "timeout_per_task": 1800.0,
        "max_retries": 1
    },
    "output": {
        "format": "json"
    }
}


@pytest.fixture
def parallel_input_file(tmp_path):
    """Create a temporary input file with parallel execution config."""
    input_file = tmp_path / "input_parallel.json"
    with open(input_file, 'w') as f:
        json.dump(BASIC_INPUT_PARALLEL, f)
    return input_file


@pytest.fixture
def sequential_input_file(tmp_path):
    """Create a temporary input file without parallel execution."""
    input_data = deepcopy(BASIC_INPUT_PARALLEL)
    input_data["execution"] = {"parallel": False}
    input_file = tmp_path / "input_sequential.json"
    with open(input_file, 'w') as f:
        json.dump(input_data, f)
    return input_file


class TestCLIParallelParsing:
    """Test parsing of parallel execution configuration from input files."""

    def test_parse_input_with_parallel_config(self, parallel_input_file):
        """Test parsing input file with parallel execution configuration."""
        from qcmanybody.cli.input_parser import parse_input_file

        cli_input = parse_input_file(parallel_input_file)

        assert cli_input.execution is not None
        assert cli_input.execution.parallel is True
        assert cli_input.execution.n_workers == 2
        assert cli_input.execution.executor_type == "multiprocessing"
        assert cli_input.execution.timeout_per_task == 1800.0
        assert cli_input.execution.max_retries == 1

    def test_parse_input_without_parallel_config(self, sequential_input_file):
        """Test parsing input file without parallel execution."""
        from qcmanybody.cli.input_parser import parse_input_file

        cli_input = parse_input_file(sequential_input_file)

        assert cli_input.execution is not None
        assert cli_input.execution.parallel is False
        # Note: When parallel=False, n_workers may be None or ignored

    def test_parse_input_with_auto_workers(self, tmp_path):
        """Test parsing input with auto-detect workers (n_workers=null)."""
        input_data = deepcopy(BASIC_INPUT_PARALLEL)
        input_data["execution"]["n_workers"] = None
        input_file = tmp_path / "input_auto.json"
        with open(input_file, 'w') as f:
            json.dump(input_data, f)

        from qcmanybody.cli.input_parser import parse_input_file

        cli_input = parse_input_file(input_file)

        assert cli_input.execution.parallel is True
        assert cli_input.execution.n_workers is None  # Should auto-detect


class TestCLIParallelExecution:
    """Test CLI execution with parallel configuration."""

    @patch('qcmanybody.cli.commands.run.format_output')
    @patch('qcmanybody.parallel.ParallelManyBodyComputer')
    @patch('qcmanybody.parallel.executors.MultiprocessingExecutor')
    def test_run_with_parallel_from_input_file(
        self,
        mock_mp_executor,
        mock_parallel_computer,
        mock_format_output,
        parallel_input_file
    ):
        """Test running calculation with parallel config from input file."""
        from qcmanybody.cli.commands.run import handle_run

        # Mock the result and formatting
        mock_result = Mock()
        mock_result.ret_energy = -5.0
        mock_parallel_computer.from_manybodyinput.return_value = mock_result
        mock_format_output.return_value = '{"ret_energy": -5.0}'

        # Create args
        args = Namespace(
            input=str(parallel_input_file),
            output=None,
            format=None,
            pretty=False,
            parallel=None,  # Not overridden
            n_workers=None  # Not overridden
        )

        # Execute
        exit_code = handle_run(args)

        # Verify success
        assert exit_code == 0

        # Verify parallel computer was called
        assert mock_parallel_computer.from_manybodyinput.called

        # Verify multiprocessing executor was created
        assert mock_mp_executor.called

    @patch('qcmanybody.cli.commands.run.format_output')
    @patch('qcmanybody.ManyBodyComputer')
    def test_run_with_sequential_from_input_file(
        self,
        mock_computer,
        mock_format_output,
        sequential_input_file
    ):
        """Test running calculation with sequential config from input file."""
        from qcmanybody.cli.commands.run import handle_run

        # Mock the result and formatting
        mock_result = Mock()
        mock_result.ret_energy = -5.0
        mock_computer.from_manybodyinput.return_value = mock_result
        mock_format_output.return_value = '{"ret_energy": -5.0}'

        # Create args
        args = Namespace(
            input=str(sequential_input_file),
            output=None,
            format=None,
            pretty=False,
            parallel=None,
            n_workers=None
        )

        # Execute
        exit_code = handle_run(args)

        # Verify success
        assert exit_code == 0

        # Verify sequential computer was called
        assert mock_computer.from_manybodyinput.called

    @patch('qcmanybody.cli.commands.run.format_output')
    @patch('qcmanybody.parallel.ParallelManyBodyComputer')
    @patch('qcmanybody.parallel.executors.MultiprocessingExecutor')
    def test_run_with_cli_parallel_override(
        self,
        mock_mp_executor,
        mock_parallel_computer,
        mock_format_output,
        sequential_input_file
    ):
        """Test that CLI arguments override input file parallel config."""
        from qcmanybody.cli.commands.run import handle_run

        # Mock the result and formatting
        mock_result = Mock()
        mock_result.ret_energy = -5.0
        mock_parallel_computer.from_manybodyinput.return_value = mock_result
        mock_format_output.return_value = '{"ret_energy": -5.0}'

        # Create args with CLI override
        args = Namespace(
            input=str(sequential_input_file),
            output=None,
            format=None,
            pretty=False,
            parallel=True,  # Override to enable parallel
            n_workers=4     # Override worker count
        )

        # Execute
        exit_code = handle_run(args)

        # Verify success
        assert exit_code == 0

        # Verify parallel computer was called (overrode sequential config)
        assert mock_parallel_computer.from_manybodyinput.called

    @patch('qcmanybody.cli.commands.run.format_output')
    @patch('qcmanybody.parallel.ParallelManyBodyComputer')
    @patch('qcmanybody.parallel.executors.SequentialExecutor')
    def test_run_with_sequential_executor_type(
        self,
        mock_seq_executor,
        mock_parallel_computer,
        mock_format_output,
        tmp_path
    ):
        """Test running with sequential executor type (for testing/debugging)."""
        # Create input with sequential executor type
        input_data = deepcopy(BASIC_INPUT_PARALLEL)
        input_data["execution"]["executor_type"] = "sequential"
        input_file = tmp_path / "input_seq_executor.json"
        with open(input_file, 'w') as f:
            json.dump(input_data, f)

        from qcmanybody.cli.commands.run import handle_run

        # Mock the result and formatting
        mock_result = Mock()
        mock_result.ret_energy = -5.0
        mock_parallel_computer.from_manybodyinput.return_value = mock_result
        mock_format_output.return_value = '{"ret_energy": -5.0}'

        # Create args
        args = Namespace(
            input=str(input_file),
            output=None,
            format=None,
            pretty=False,
            parallel=None,
            n_workers=None
        )

        # Execute
        exit_code = handle_run(args)

        # Verify success
        assert exit_code == 0

        # Verify sequential executor was used
        assert mock_seq_executor.called


class TestCLIMainArguments:
    """Test CLI main.py argument parsing for parallel execution."""

    def test_parallel_arguments_defined(self):
        """Test that parallel execution arguments are defined in CLI."""
        from qcmanybody.cli.main import create_parser

        parser = create_parser()

        # Parse with parallel arguments
        args = parser.parse_args(['run', 'input.json', '--parallel', '--n-workers', '4'])

        assert args.command == 'run'
        assert args.parallel is True
        assert args.n_workers == 4

    def test_no_parallel_flag(self):
        """Test --no-parallel flag disables parallel execution."""
        from qcmanybody.cli.main import create_parser

        parser = create_parser()

        # Parse with --no-parallel
        args = parser.parse_args(['run', 'input.json', '--no-parallel'])

        assert args.parallel is False

    def test_default_parallel_is_none(self):
        """Test default parallel value is None (use input file config)."""
        from qcmanybody.cli.main import create_parser

        parser = create_parser()

        # Parse without parallel arguments
        args = parser.parse_args(['run', 'input.json'])

        # Default should be None to allow input file to control
        assert args.parallel is None

    def test_n_workers_validation(self):
        """Test n_workers argument accepts positive integers."""
        from qcmanybody.cli.main import create_parser

        parser = create_parser()

        # Valid worker counts
        for n in [1, 2, 4, 8, 16]:
            args = parser.parse_args(['run', 'input.json', '--n-workers', str(n)])
            assert args.n_workers == n


class TestExecutorConfigCreation:
    """Test creation of ExecutorConfig from CLI input."""

    def test_executor_config_from_cli_input(self, parallel_input_file):
        """Test that ExecutorConfig is properly created from CLI input."""
        from qcmanybody.cli.input_parser import parse_input_file
        from qcmanybody.parallel import ExecutorConfig

        cli_input = parse_input_file(parallel_input_file)
        exec_config = cli_input.execution

        # Verify the execution config has expected values
        assert exec_config.n_workers == 2
        assert exec_config.timeout_per_task == 1800.0
        assert exec_config.max_retries == 1

        # Create ExecutorConfig and verify it works
        config = ExecutorConfig(
            n_workers=exec_config.n_workers,
            timeout_per_task=exec_config.timeout_per_task,
            max_retries=exec_config.max_retries
        )

        assert config.n_workers == 2
        assert config.timeout_per_task == 1800.0
        assert config.max_retries == 1

    def test_executor_config_with_defaults(self):
        """Test ExecutorConfig with CLI defaults."""
        from qcmanybody.parallel import ExecutorConfig

        config = ExecutorConfig(
            n_workers=None,  # Auto-detect
            timeout_per_task=3600.0,  # Default 1 hour
            max_retries=2  # Default
        )

        assert config.n_workers is None
        assert config.timeout_per_task == 3600.0
        assert config.max_retries == 2


class TestCLIErrorHandling:
    """Test error handling in CLI with parallel execution."""

    def test_missing_qcengine_error(self, parallel_input_file):
        """Test graceful error when qcengine is not installed."""
        from qcmanybody.cli.commands.run import handle_run

        # Mock ImportError for qcengine
        with patch('qcmanybody.parallel.ParallelManyBodyComputer') as mock_computer:
            mock_computer.from_manybodyinput.side_effect = ImportError("No module named 'qcengine'")

            args = Namespace(
                input=str(parallel_input_file),
                output=None,
                format=None,
                pretty=False,
                parallel=None,
                n_workers=None
            )

            exit_code = handle_run(args)

            # Should return non-zero exit code
            assert exit_code != 0

    def test_invalid_input_file(self, tmp_path):
        """Test error handling for invalid input file."""
        from qcmanybody.cli.commands.run import handle_run

        # Create invalid input file
        invalid_file = tmp_path / "invalid.json"
        with open(invalid_file, 'w') as f:
            f.write("{ invalid json ")

        args = Namespace(
            input=str(invalid_file),
            output=None,
            format=None,
            pretty=False,
            parallel=None,
            n_workers=None
        )

        exit_code = handle_run(args)

        # Should return non-zero exit code
        assert exit_code != 0

    def test_calculation_failure_returns_error(self, parallel_input_file):
        """Test that calculation failures return non-zero exit code."""
        from qcmanybody.cli.commands.run import handle_run

        # Mock calculation failure
        with patch('qcmanybody.parallel.ParallelManyBodyComputer') as mock_computer:
            mock_computer.from_manybodyinput.side_effect = RuntimeError("Calculation failed")

            args = Namespace(
                input=str(parallel_input_file),
                output=None,
                format=None,
                pretty=False,
                parallel=None,
                n_workers=None
            )

            exit_code = handle_run(args)

            # Should return non-zero exit code
            assert exit_code != 0


@pytest.mark.integration
class TestEndToEndCLI:
    """End-to-end tests for CLI with parallel execution (require qcengine)."""

    @pytest.mark.skip(reason="Requires qcengine with QC program installed")
    def test_full_parallel_execution(self, parallel_input_file):
        """Test full parallel execution with real QC calculation."""
        # This test would require actual qcengine + psi4/other QC program
        # Skip by default, enable for full integration testing
        pass

    @pytest.mark.skip(reason="Requires qcengine with QC program installed")
    def test_parallel_vs_sequential_results_match(self, tmp_path):
        """Test that parallel and sequential execution produce identical results."""
        # This test would verify that parallel execution produces same results
        # Skip by default, enable for full integration testing
        pass
