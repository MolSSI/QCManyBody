"""
Performance benchmarks for CLI.

These tests measure the performance of CLI operations to ensure they're efficient.
"""

import json
import subprocess
import tempfile
import time
from pathlib import Path

import pytest


def create_test_input():
    """Create a test input file for benchmarking."""
    return {
        "schema_name": "qcmanybody_cli_input",
        "schema_version": 1,
        "molecule": {
            "source": "inline",
            "inline": {
                "symbols": ["He", "He", "He"],
                "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0], [0.0, 0.0, 6.0]],
                "fragments": [[0], [1], [2]],
                "units": "angstrom",
            },
        },
        "calculation": {
            "type": "single",
            "single": {"driver": "energy", "method": "hf", "basis": "sto-3g", "program": "psi4"},
        },
        "bsse": {"type": ["cp"]},
        "manybody": {"max_nbody": 3},
    }


def test_cli_validation_performance():
    """Benchmark validation command performance."""
    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        # Warm up
        subprocess.run(["qcmanybody", "validate", temp_path], capture_output=True)

        # Benchmark
        iterations = 10
        start = time.time()
        for _ in range(iterations):
            result = subprocess.run(["qcmanybody", "validate", temp_path], capture_output=True)
            assert result.returncode == 0
        end = time.time()

        avg_time = (end - start) / iterations
        print(f"\nValidation average time: {avg_time:.3f}s ({iterations} iterations)")

        # Validation should be fast (< 2 seconds per call)
        assert avg_time < 2.0, f"Validation too slow: {avg_time:.3f}s"

    finally:
        Path(temp_path).unlink()


def test_cli_plan_performance():
    """Benchmark plan command performance."""
    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        # Warm up
        subprocess.run(["qcmanybody", "plan", temp_path], capture_output=True)

        # Benchmark
        iterations = 10
        start = time.time()
        for _ in range(iterations):
            result = subprocess.run(["qcmanybody", "plan", temp_path], capture_output=True)
            assert result.returncode == 0
        end = time.time()

        avg_time = (end - start) / iterations
        print(f"\nPlan average time: {avg_time:.3f}s ({iterations} iterations)")

        # Planning should be fast (< 2 seconds per call)
        assert avg_time < 2.0, f"Planning too slow: {avg_time:.3f}s"

    finally:
        Path(temp_path).unlink()


def test_cli_convert_performance():
    """Benchmark convert command performance."""
    pytest.importorskip("yaml")

    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f_in:
        json.dump(json_data, f_in)
        input_path = f_in.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f_out:
        output_path = f_out.name

    try:
        # Warm up
        subprocess.run(["qcmanybody", "convert", input_path, output_path], capture_output=True)

        # Benchmark
        iterations = 10
        start = time.time()
        for _ in range(iterations):
            result = subprocess.run(["qcmanybody", "convert", input_path, output_path], capture_output=True)
            assert result.returncode == 0
        end = time.time()

        avg_time = (end - start) / iterations
        print(f"\nConvert average time: {avg_time:.3f}s ({iterations} iterations)")

        # Conversion should be fast (< 2 seconds per call)
        assert avg_time < 2.0, f"Conversion too slow: {avg_time:.3f}s"

    finally:
        Path(input_path).unlink()
        Path(output_path).unlink(missing_ok=True)


def test_cli_startup_overhead():
    """Measure CLI startup overhead."""
    # Test --version command (minimal work)
    iterations = 20
    start = time.time()
    for _ in range(iterations):
        subprocess.run(["qcmanybody", "--version"], capture_output=True)
    end = time.time()

    avg_time = (end - start) / iterations
    print(f"\nCLI startup overhead: {avg_time:.3f}s ({iterations} iterations)")

    # Startup should be fast (< 1 second)
    assert avg_time < 1.0, f"CLI startup too slow: {avg_time:.3f}s"


def test_python_api_vs_cli_overhead():
    """Compare Python API vs CLI for input parsing."""
    from qcmanybody.cli.converter import convert_to_manybody_input
    from qcmanybody.cli.input_parser import parse_input_file

    json_data = create_test_input()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        # Benchmark Python API
        iterations = 50
        start = time.time()
        for _ in range(iterations):
            cli_input = parse_input_file(temp_path)
            mb_input = convert_to_manybody_input(cli_input, temp_path)
        end = time.time()
        python_time = (end - start) / iterations

        # Benchmark CLI (validation only, no execution)
        start = time.time()
        for _ in range(iterations):
            subprocess.run(["qcmanybody", "validate", temp_path], capture_output=True)
        end = time.time()
        cli_time = (end - start) / iterations

        overhead = cli_time - python_time
        overhead_percent = (overhead / python_time) * 100

        print(f"\nPython API parse+convert: {python_time:.4f}s")
        print(f"CLI validate: {cli_time:.4f}s")
        print(f"CLI overhead: {overhead:.4f}s ({overhead_percent:.1f}%)")

        # CLI overhead should be reasonable (< 500ms more than direct Python)
        assert overhead < 0.5, f"CLI overhead too high: {overhead:.3f}s"

    finally:
        Path(temp_path).unlink()


def test_input_file_size_scaling():
    """Test performance with different input file sizes."""
    results = []

    # Test with different numbers of fragments (scaling complexity)
    for n_fragments in [2, 3, 5]:
        json_data = {
            "schema_name": "qcmanybody_cli_input",
            "schema_version": 1,
            "molecule": {
                "source": "inline",
                "inline": {
                    "symbols": ["He"] * n_fragments,
                    "geometry": [[0.0, 0.0, float(i * 3)] for i in range(n_fragments)],
                    "fragments": [[i] for i in range(n_fragments)],
                    "units": "angstrom",
                },
            },
            "calculation": {
                "type": "single",
                "single": {"driver": "energy", "method": "hf", "basis": "sto-3g", "program": "psi4"},
            },
            "bsse": {"type": ["cp"]},
            "manybody": {"max_nbody": min(n_fragments, 3)},
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(json_data, f)
            temp_path = f.name

        try:
            iterations = 10
            start = time.time()
            for _ in range(iterations):
                subprocess.run(["qcmanybody", "plan", temp_path], capture_output=True)
            end = time.time()

            avg_time = (end - start) / iterations
            results.append((n_fragments, avg_time))
            print(f"\n{n_fragments} fragments: {avg_time:.3f}s")

        finally:
            Path(temp_path).unlink()

    # Performance should scale reasonably (not exponentially)
    # Allow up to 50% increase per additional fragment for planning
    for i in range(1, len(results)):
        ratio = results[i][1] / results[i - 1][1]
        assert ratio < 1.5, f"Performance scaling too steep: {results[i - 1][0]}->{results[i][0]} fragments: {ratio:.2f}x"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
