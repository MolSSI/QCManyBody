"""
Tests for CLI example input files.

Validates that all example input files in examples/cli/ are valid and can be parsed.
"""

import subprocess
from pathlib import Path

import pytest


def find_example_files():
    """Find all example input files in examples/cli/ directory."""
    examples_dir = Path(__file__).parent.parent.parent / "examples" / "cli"
    if not examples_dir.exists():
        return []

    # Find all JSON and YAML files
    json_files = list(examples_dir.glob("*.json"))
    yaml_files = list(examples_dir.glob("*.yaml"))
    yml_files = list(examples_dir.glob("*.yml"))

    all_files = json_files + yaml_files + yml_files

    # Filter out non-input files
    input_files = [f for f in all_files if not f.stem.endswith(("_output", "_result"))]

    return input_files


@pytest.mark.parametrize("example_file", find_example_files(), ids=lambda f: f.name)
def test_example_file_validates(example_file):
    """Test that each example file validates successfully."""
    result = subprocess.run(
        ["qcmanybody", "validate", str(example_file)],
        capture_output=True,
        text=True,
    )

    # Check that validation succeeds
    assert result.returncode == 0, (
        f"Validation failed for {example_file.name}\n"
        f"stdout: {result.stdout}\n"
        f"stderr: {result.stderr}"
    )
    assert "Validation successful" in result.stdout or "✓" in result.stdout


@pytest.mark.parametrize("example_file", find_example_files(), ids=lambda f: f.name)
def test_example_file_plan(example_file):
    """Test that each example file can generate an execution plan."""
    result = subprocess.run(
        ["qcmanybody", "plan", str(example_file)],
        capture_output=True,
        text=True,
    )

    # Check that plan generation succeeds
    assert result.returncode == 0, (
        f"Plan generation failed for {example_file.name}\n"
        f"stdout: {result.stdout}\n"
        f"stderr: {result.stderr}"
    )
    # Should show execution plan output
    assert ("Execution Plan" in result.stdout or
            "Task" in result.stdout or
            "Molecular System" in result.stdout)


def test_examples_directory_exists():
    """Test that examples directory exists and has files."""
    examples_dir = Path(__file__).parent.parent.parent / "examples" / "cli"
    assert examples_dir.exists(), f"Examples directory not found at {examples_dir}"

    example_files = find_example_files()
    assert len(example_files) > 0, "No example files found in examples/cli/"

    print(f"\nFound {len(example_files)} example files:")
    for f in example_files:
        print(f"  - {f.name}")


def test_basic_energy_example():
    """Specific test for the basic energy example."""
    examples_dir = Path(__file__).parent.parent.parent / "examples" / "cli"
    basic_energy = examples_dir / "01_basic_energy.json"

    if not basic_energy.exists():
        pytest.skip("01_basic_energy.json not found")

    # Test validation
    result = subprocess.run(
        ["qcmanybody", "validate", str(basic_energy)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0

    # Test plan
    result = subprocess.run(
        ["qcmanybody", "plan", str(basic_energy)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0


def test_multilevel_example():
    """Specific test for the multilevel example."""
    examples_dir = Path(__file__).parent.parent.parent / "examples" / "cli"
    multilevel = examples_dir / "03_multilevel.json"

    if not multilevel.exists():
        pytest.skip("03_multilevel.json not found")

    # Test validation
    result = subprocess.run(
        ["qcmanybody", "validate", str(multilevel)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0

    # Test plan
    result = subprocess.run(
        ["qcmanybody", "plan", str(multilevel)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    # Multi-level should show different methods
    assert "ccsd" in result.stdout.lower() or "mp2" in result.stdout.lower()
