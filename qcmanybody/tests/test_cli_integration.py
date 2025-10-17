"""
Integration tests for CLI commands.

These tests verify end-to-end functionality of CLI commands without running actual QC calculations.
"""

import json
import subprocess
import tempfile
from pathlib import Path
from unittest import mock

import pytest


def create_test_input(driver="energy", max_nbody=2):
    """Create a basic test input file."""
    return {
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
                "units": "angstrom",
            },
        },
        "calculation": {
            "type": "single",
            "single": {"driver": driver, "method": "hf", "basis": "sto-3g", "program": "psi4"},
        },
        "manybody": {
            "max_nbody": max_nbody,
        },
    }


def test_cli_validate_valid_input():
    """Test validate command with valid input."""
    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        result = subprocess.run(
            ["qcmanybody", "validate", temp_path],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Validation successful" in result.stdout or "✓" in result.stdout
    finally:
        Path(temp_path).unlink()


def test_cli_validate_invalid_input():
    """Test validate command with invalid input."""
    json_data = {
        "schema_name": "qcmanybody_cli_input",
        "schema_version": 1,
        # Missing required fields
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        result = subprocess.run(
            ["qcmanybody", "validate", temp_path],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "Validation failed" in result.stdout or "Error" in result.stdout or "error" in result.stderr
    finally:
        Path(temp_path).unlink()


def test_cli_validate_show_schema():
    """Test validate command with --show-schema option."""
    result = subprocess.run(
        ["qcmanybody", "validate", "--show-schema"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "qcmanybody_cli_input" in result.stdout
    assert "schema" in result.stdout.lower()


def test_cli_plan():
    """Test plan command."""
    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        result = subprocess.run(
            ["qcmanybody", "plan", temp_path],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Execution Plan" in result.stdout or "Task" in result.stdout
        # Should show molecular system info
        assert "Molecular System" in result.stdout or "fragments" in result.stdout.lower()
    finally:
        Path(temp_path).unlink()


def test_cli_plan_show_tasks():
    """Test plan command with --show-tasks option."""
    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        result = subprocess.run(
            ["qcmanybody", "plan", temp_path, "--show-tasks"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        # Should show detailed task information
        assert "body" in result.stdout.lower() or "calculation" in result.stdout.lower()
    finally:
        Path(temp_path).unlink()


def test_cli_convert_json_to_yaml():
    """Test convert command from JSON to YAML."""
    pytest.importorskip("yaml")  # Skip if PyYAML not available

    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f_in:
        json.dump(json_data, f_in)
        input_path = f_in.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f_out:
        output_path = f_out.name

    try:
        result = subprocess.run(
            ["qcmanybody", "convert", input_path, output_path],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert Path(output_path).exists()

        # Verify YAML content
        import yaml
        with open(output_path, "r") as f:
            yaml_data = yaml.safe_load(f)
        assert yaml_data["schema_name"] == "qcmanybody_cli_input"
        assert yaml_data["calculation"]["type"] == "single"
    finally:
        Path(input_path).unlink()
        Path(output_path).unlink(missing_ok=True)


def test_cli_convert_yaml_to_json():
    """Test convert command from YAML to JSON."""
    pytest.importorskip("yaml")  # Skip if PyYAML not available

    import yaml

    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f_in:
        yaml.dump(json_data, f_in)
        input_path = f_in.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f_out:
        output_path = f_out.name

    try:
        result = subprocess.run(
            ["qcmanybody", "convert", input_path, output_path],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert Path(output_path).exists()

        # Verify JSON content
        with open(output_path, "r") as f:
            json_data_out = json.load(f)
        assert json_data_out["schema_name"] == "qcmanybody_cli_input"
        assert json_data_out["calculation"]["type"] == "single"
    finally:
        Path(input_path).unlink()
        Path(output_path).unlink(missing_ok=True)


def test_cli_main_help():
    """Test main help output."""
    result = subprocess.run(
        ["qcmanybody", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "run" in result.stdout
    assert "plan" in result.stdout
    assert "validate" in result.stdout
    assert "convert" in result.stdout


def test_cli_run_help():
    """Test run command help output."""
    result = subprocess.run(
        ["qcmanybody", "run", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "input" in result.stdout.lower()
    assert "output" in result.stdout.lower()


def test_cli_version():
    """Test version output."""
    result = subprocess.run(
        ["qcmanybody", "--version"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    # Should print version info
    assert len(result.stdout.strip()) > 0 or len(result.stderr.strip()) > 0


def test_cli_run_returns_manybody_result():
    """Test that run command properly handles ManyBodyResult from ManyBodyComputer.from_manybodyinput."""
    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f_out:
        output_path = f_out.name

    try:
        # Create a mock ManyBodyResult
        from qcmanybody.models.v1 import ManyBodyResult, ManyBodyInput

        mock_result = mock.Mock(spec=ManyBodyResult)
        mock_result.dict.return_value = {
            "input_data": {"schema_name": "qcschema_manybodyinput"},
            "properties": {
                "return_energy": -5.0,
                "calcinfo_nmbe": 3,
            },
            "success": True,
        }

        # Mock the ManyBodyComputer at the source module level
        with mock.patch("qcmanybody.ManyBodyComputer") as mock_computer:
            mock_computer.from_manybodyinput.return_value = mock_result

            # Call run command via main.main to test in-process
            from qcmanybody.cli.main import main

            result = main(["run", temp_path, "-o", output_path])

            # Should succeed
            assert result == 0

            # Verify output file was written
            assert Path(output_path).exists()

            # Verify output contains the mocked result
            with open(output_path, "r") as f:
                output_data = json.load(f)

            assert "results" in output_data
            assert output_data["results"]["success"] is True
            assert output_data["results"]["properties"]["return_energy"] == -5.0

    finally:
        Path(temp_path).unlink()
        Path(output_path).unlink(missing_ok=True)


def test_cli_run_format_precedence():
    """Test that CLI format flag overrides schema format."""
    json_data = create_test_input()
    json_data["output"] = {"format": "text"}  # Schema says text

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f_out:
        output_path = f_out.name

    try:
        mock_result = mock.Mock()
        mock_result.dict.return_value = {"success": True, "properties": {"return_energy": -5.0}}

        with mock.patch("qcmanybody.ManyBodyComputer") as mock_computer:
            mock_computer.from_manybodyinput.return_value = mock_result

            from qcmanybody.cli.main import main

            # CLI specifies json, should override schema's text
            result = main(["run", temp_path, "-o", output_path, "--format", "json"])

            assert result == 0
            with open(output_path, "r") as f:
                output_data = json.load(f)  # Should be valid JSON, not text
            assert "results" in output_data

    finally:
        Path(temp_path).unlink()
        Path(output_path).unlink(missing_ok=True)


def test_cli_run_output_path_precedence():
    """Test that CLI output path overrides schema output path."""
    json_data = create_test_input()
    json_data["output"] = {"file": "schema_output.json"}  # Schema says one path

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f_out:
        cli_output_path = f_out.name

    try:
        mock_result = mock.Mock()
        mock_result.dict.return_value = {"success": True, "properties": {"return_energy": -5.0}}

        with mock.patch("qcmanybody.ManyBodyComputer") as mock_computer:
            mock_computer.from_manybodyinput.return_value = mock_result

            from qcmanybody.cli.main import main

            # CLI specifies different path, should override schema
            result = main(["run", temp_path, "-o", cli_output_path])

            assert result == 0
            # Should write to CLI path, not schema path
            assert Path(cli_output_path).exists()
            assert not Path("schema_output.json").exists()

    finally:
        Path(temp_path).unlink()
        Path(cli_output_path).unlink(missing_ok=True)
        Path("schema_output.json").unlink(missing_ok=True)


def test_cli_run_schema_only_defaults():
    """Test that schema-only settings are used when CLI flags not specified."""
    json_data = create_test_input()
    json_data["output"] = {"file": "output_from_schema.json", "format": "json"}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        mock_result = mock.Mock()
        mock_result.dict.return_value = {"success": True, "properties": {"return_energy": -5.0}}

        with mock.patch("qcmanybody.ManyBodyComputer") as mock_computer:
            mock_computer.from_manybodyinput.return_value = mock_result

            from qcmanybody.cli.main import main

            # No CLI flags, should use schema settings
            result = main(["run", temp_path])

            assert result == 0
            # Should write to schema path
            assert Path("output_from_schema.json").exists()

    finally:
        Path(temp_path).unlink()
        Path("output_from_schema.json").unlink(missing_ok=True)


# Note: Full end-to-end tests with actual QC calculations require QC programs (like Psi4)
# to be installed. Those would be tested separately with appropriate test markers.
