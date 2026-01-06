"""
Unit tests for CLI converter module.
"""

import json
import tempfile
from pathlib import Path

import pytest

from qcmanybody.cli.converter import convert_to_manybody_input
from qcmanybody.cli.input_parser import parse_input_file


def test_convert_single_level():
    """Test converting a single-level calculation."""
    json_data = {
        "schema_name": "qcmanybody_cli_input",
        "schema_version": 1,
        "molecule": {
            "source": "inline",
            "inline": {
                "symbols": ["He", "He"],
                "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]],
                "fragments": [[0], [1]],
                "molecular_charge": 0.0,
                "molecular_multiplicity": 1,
                "units": "angstrom",
            },
        },
        "calculation": {
            "type": "single",
            "single": {"driver": "energy", "method": "hf", "basis": "sto-3g", "program": "psi4"},
        },
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        cli_input = parse_input_file(temp_path)
        mb_input = convert_to_manybody_input(cli_input, input_file_path=temp_path)

        # Check molecule
        assert len(mb_input.molecule.symbols) == 2
        assert len(mb_input.molecule.fragments) == 2
        assert list(mb_input.molecule.symbols) == ["He", "He"]

        # Check specification
        assert mb_input.specification.driver.value == "energy"
        assert "(auto)" in mb_input.specification.specification
        atomic_spec = mb_input.specification.specification["(auto)"]
        assert atomic_spec.model.method == "hf"
        assert atomic_spec.model.basis == "sto-3g"
        assert atomic_spec.program == "psi4"

    finally:
        Path(temp_path).unlink()


def test_convert_multilevel():
    """Test converting a multi-level calculation."""
    json_data = {
        "schema_name": "qcmanybody_cli_input",
        "schema_version": 1,
        "molecule": {
            "source": "inline",
            "inline": {
                "symbols": ["He", "He", "He"],
                "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 4.0]],
                "fragments": [[0], [1], [2]],
                "molecular_charge": 0.0,
                "molecular_multiplicity": 1,
                "units": "angstrom",
            },
        },
        "calculation": {
            "type": "multi",
            "multi": {
                "driver": "energy",
                "levels": {
                    "1": {"method": "ccsd", "basis": "cc-pvdz", "program": "psi4"},
                    "2": {"method": "mp2", "basis": "cc-pvdz", "program": "psi4"},
                    "3": {"method": "hf", "basis": "sto-3g", "program": "psi4"},
                },
            },
        },
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        cli_input = parse_input_file(temp_path)
        mb_input = convert_to_manybody_input(cli_input, input_file_path=temp_path)

        # Check that multi-level specs are created with method/basis keys
        assert len(mb_input.specification.specification) == 3
        assert "ccsd/cc-pvdz" in mb_input.specification.specification
        assert "mp2/cc-pvdz" in mb_input.specification.specification
        assert "hf/sto-3g" in mb_input.specification.specification

        # Check level 1 specification
        spec1 = mb_input.specification.specification["ccsd/cc-pvdz"]
        assert spec1.model.method == "ccsd"
        assert spec1.model.basis == "cc-pvdz"

        # Check levels mapping in keywords
        assert mb_input.specification.keywords.levels is not None
        assert mb_input.specification.keywords.levels[1] == "ccsd/cc-pvdz"
        assert mb_input.specification.keywords.levels[2] == "mp2/cc-pvdz"
        assert mb_input.specification.keywords.levels[3] == "hf/sto-3g"

    finally:
        Path(temp_path).unlink()


def test_convert_bsse_types():
    """Test that BSSE types are correctly mapped."""
    json_data = {
        "schema_name": "qcmanybody_cli_input",
        "schema_version": 1,
        "molecule": {
            "source": "inline",
            "inline": {
                "symbols": ["He", "He"],
                "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]],
                "fragments": [[0], [1]],
                "molecular_charge": 0.0,
                "molecular_multiplicity": 1,
                "units": "angstrom",
            },
        },
        "calculation": {
            "type": "single",
            "single": {"driver": "energy", "method": "hf", "basis": "sto-3g", "program": "psi4"},
        },
        "bsse": {"type": ["cp", "vmfc"]},
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        cli_input = parse_input_file(temp_path)
        mb_input = convert_to_manybody_input(cli_input, input_file_path=temp_path)

        # Check BSSE types
        bsse_types = mb_input.specification.keywords.bsse_type
        assert len(bsse_types) == 2
        bsse_values = [bt.value for bt in bsse_types]
        assert "cp" in bsse_values
        assert "vmfc" in bsse_values

    finally:
        Path(temp_path).unlink()


def test_convert_manybody_keywords():
    """Test that many-body keywords are correctly set."""
    json_data = {
        "schema_name": "qcmanybody_cli_input",
        "schema_version": 1,
        "molecule": {
            "source": "inline",
            "inline": {
                "symbols": ["He", "He", "He"],
                "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 4.0]],
                "fragments": [[0], [1], [2]],
                "molecular_charge": 0.0,
                "molecular_multiplicity": 1,
                "units": "angstrom",
            },
        },
        "calculation": {
            "type": "single",
            "single": {"driver": "energy", "method": "hf", "basis": "sto-3g", "program": "psi4"},
        },
        "manybody": {
            "max_nbody": 2,
            "return_total_data": True,
            "supersystem_ie_only": False,
        },
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        cli_input = parse_input_file(temp_path)
        mb_input = convert_to_manybody_input(cli_input, input_file_path=temp_path)

        # Check many-body keywords
        assert mb_input.specification.keywords.max_nbody == 2
        assert mb_input.specification.keywords.return_total_data == True
        assert mb_input.specification.keywords.supersystem_ie_only == False

    finally:
        Path(temp_path).unlink()
