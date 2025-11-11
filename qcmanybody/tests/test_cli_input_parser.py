"""
Unit tests for CLI input parser module.
"""

import json
import tempfile
from pathlib import Path

import pytest

from qcmanybody.cli.input_parser import parse_input_file


def test_parse_json_basic():
    """Test parsing a basic JSON input file."""
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
            "single": {
                "driver": "energy",
                "method": "hf",
                "basis": "sto-3g",
                "program": "psi4",
            },
        },
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        cli_input = parse_input_file(temp_path)
        assert cli_input.schema_version == 1
        assert cli_input.molecule.source.value == "inline"
        assert len(cli_input.molecule.inline.symbols) == 2
        assert cli_input.calculation.single.driver.value == "energy"
        assert cli_input.calculation.single.method == "hf"
    finally:
        Path(temp_path).unlink()


def test_parse_yaml_basic():
    """Test parsing a basic YAML input file."""
    pytest.importorskip("yaml")

    yaml_content = """
schema_name: qcmanybody_cli_input
schema_version: 1
molecule:
  source: inline
  inline:
    symbols: [He, He]
    geometry: [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]]
    fragments: [[0], [1]]
    molecular_charge: 0.0
    molecular_multiplicity: 1
    units: angstrom
calculation:
  type: single
  single:
    driver: energy
    method: hf
    basis: sto-3g
    program: psi4
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        temp_path = f.name

    try:
        cli_input = parse_input_file(temp_path)
        assert cli_input.schema_version == 1
        assert cli_input.molecule.source.value == "inline"
        assert len(cli_input.molecule.inline.symbols) == 2
    finally:
        Path(temp_path).unlink()


def test_parse_invalid_json():
    """Test that invalid JSON raises an error."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        f.write("{invalid json")
        temp_path = f.name

    try:
        with pytest.raises(Exception):
            parse_input_file(temp_path)
    finally:
        Path(temp_path).unlink()


def test_parse_missing_required_field():
    """Test that missing required fields raise validation errors."""
    json_data = {
        "schema_name": "qcmanybody_cli_input",
        "schema_version": 1,
        # Missing molecule and calculation
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(json_data, f)
        temp_path = f.name

    try:
        with pytest.raises(Exception):
            parse_input_file(temp_path)
    finally:
        Path(temp_path).unlink()


def test_parse_multilevel_calculation():
    """Test parsing a multi-level calculation."""
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
        assert cli_input.calculation.type == "multi"
        assert cli_input.calculation.multi is not None
        assert len(cli_input.calculation.multi.levels) == 3
        # Levels are stored as integer keys after validation
        assert 1 in cli_input.calculation.multi.levels
        assert cli_input.calculation.multi.levels[1]["method"] == "ccsd"
        assert cli_input.calculation.multi.levels[1]["basis"] == "cc-pvdz"
        assert cli_input.calculation.multi.levels[1]["program"] == "psi4"
    finally:
        Path(temp_path).unlink()


def test_parse_bsse_types():
    """Test parsing different BSSE type specifications."""
    # Single BSSE type
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
        assert len(cli_input.bsse.type) == 2
        assert cli_input.bsse.type[0].value == "cp"
        assert cli_input.bsse.type[1].value == "vmfc"
    finally:
        Path(temp_path).unlink()


def test_parse_file_not_found():
    """Test that non-existent file raises appropriate error."""
    with pytest.raises(Exception):
        parse_input_file("/nonexistent/file.json")


def test_molecule_units_angstrom():
    """Test that angstrom units are properly preserved in inline molecules."""
    from qcmanybody.cli.molecule_loader import create_inline_molecule
    from qcmanybody.cli.schemas.input_schema import InlineMoleculeSchema, MoleculeSchema, MoleculeSourceEnum

    inline_schema = InlineMoleculeSchema(
        symbols=["He", "He"],
        geometry=[[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]],
        fragments=[[0], [1]],
        units="angstrom",
    )
    mol_schema = MoleculeSchema(source=MoleculeSourceEnum.inline, inline=inline_schema)

    mol = create_inline_molecule(mol_schema)

    # Verify molecule was created with angstrom units
    # QCElemental stores geometry internally in Bohr, but we can check the distance
    # Distance should be 2.0 Angstroms = ~3.78 Bohr
    import numpy as np

    dist_bohr = np.linalg.norm(np.array(mol.geometry[1]) - np.array(mol.geometry[0]))
    expected_dist_bohr = 2.0 * 1.8897259886  # 2.0 Angstrom to Bohr conversion
    assert abs(dist_bohr - expected_dist_bohr) < 0.001, f"Distance {dist_bohr} != {expected_dist_bohr}"


def test_molecule_units_bohr():
    """Test that bohr units are properly preserved in inline molecules."""
    from qcmanybody.cli.molecule_loader import create_inline_molecule
    from qcmanybody.cli.schemas.input_schema import InlineMoleculeSchema, MoleculeSchema, MoleculeSourceEnum

    inline_schema = InlineMoleculeSchema(
        symbols=["He", "He"],
        geometry=[[0.0, 0.0, 0.0], [0.0, 0.0, 3.78]],  # ~2 Angstroms in Bohr
        fragments=[[0], [1]],
        units="bohr",
    )
    mol_schema = MoleculeSchema(source=MoleculeSourceEnum.inline, inline=inline_schema)

    mol = create_inline_molecule(mol_schema)

    # Verify molecule was created with bohr units
    # QCElemental stores geometry internally in Bohr
    import numpy as np

    dist_bohr = np.linalg.norm(np.array(mol.geometry[1]) - np.array(mol.geometry[0]))
    expected_dist_bohr = 3.78
    assert abs(dist_bohr - expected_dist_bohr) < 0.001, f"Distance {dist_bohr} != {expected_dist_bohr}"


def test_xyz_file_units():
    """Test that XYZ files are properly loaded with angstrom units."""
    from qcmanybody.cli.molecule_loader import load_xyz_file

    # Create temporary XYZ file
    xyz_content = """2
Water dimer fragment
He 0.0 0.0 0.0
He 0.0 0.0 2.0
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
        f.write(xyz_content)
        temp_path = f.name

    try:
        mol = load_xyz_file(temp_path, fragments=[[0], [1]])

        # Verify distance is 2.0 Angstroms = ~3.78 Bohr in internal representation
        import numpy as np

        dist_bohr = np.linalg.norm(np.array(mol.geometry[1]) - np.array(mol.geometry[0]))
        expected_dist_bohr = 2.0 * 1.8897259886
        assert abs(dist_bohr - expected_dist_bohr) < 0.001, f"Distance {dist_bohr} != {expected_dist_bohr}"
    finally:
        Path(temp_path).unlink()
