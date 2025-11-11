"""
QCManyBody CLI Input File Schema

Defines Pydantic models for validating CLI input files (JSON/YAML).
This provides a user-friendly interface that is then converted to the internal ManyBodyInput format.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

from pydantic.v1 import BaseModel, Field, validator

# ==================== Enumerations ====================


class DriverEnum(str, Enum):
    """Type of quantum chemistry calculation."""

    energy = "energy"
    gradient = "gradient"
    hessian = "hessian"


class BsseTypeEnum(str, Enum):
    """BSSE correction methods."""

    nocp = "nocp"  # No counterpoise correction
    cp = "cp"  # Counterpoise corrected
    vmfc = "vmfc"  # Valiron-Mayer function counterpoise
    # Aliases
    ssfc = "cp"
    mbe = "nocp"
    none = "nocp"


class MoleculeSourceEnum(str, Enum):
    """How the molecule is specified."""

    inline = "inline"  # Directly in input file
    xyz = "xyz"  # From XYZ file
    pdb = "pdb"  # From PDB file
    qcschema = "qcschema"  # From QCSchema JSON file


class OutputFormatEnum(str, Enum):
    """Output file formats."""

    json = "json"
    yaml = "yaml"
    text = "text"


# ==================== Molecule Specification ====================


class InlineMoleculeSchema(BaseModel):
    """Inline molecule specification."""

    symbols: List[str] = Field(
        ...,
        description="List of atomic symbols (e.g., ['He', 'He', 'He']).",
    )
    geometry: List[List[float]] = Field(
        ...,
        description="Atomic coordinates in Angstroms or Bohr. List of [x, y, z] for each atom.",
    )
    fragments: List[List[int]] = Field(
        ...,
        description="Fragment definitions as lists of 0-based atom indices (e.g., [[0], [1], [2]]).",
    )
    molecular_charge: Optional[float] = Field(
        0.0,
        description="Total molecular charge.",
    )
    molecular_multiplicity: Optional[int] = Field(
        1,
        description="Spin multiplicity (2S+1).",
    )
    units: Optional[Literal["angstrom", "bohr"]] = Field(
        "angstrom",
        description="Units for geometry specification.",
    )
    fragment_charges: Optional[List[float]] = Field(
        None,
        description="Charge of each fragment. If not specified, assumed to be 0 for all fragments.",
    )
    fragment_multiplicities: Optional[List[int]] = Field(
        None,
        description="Multiplicity of each fragment. If not specified, assumed to be 1 for all fragments.",
    )

    @validator("geometry")
    @classmethod
    def validate_geometry(cls, v, values):
        """Ensure geometry has correct shape."""
        if "symbols" in values and len(v) != len(values["symbols"]):
            raise ValueError(f"Geometry has {len(v)} coordinates but {len(values['symbols'])} atoms specified")
        for coord in v:
            if len(coord) != 3:
                raise ValueError(f"Each coordinate must have exactly 3 values (x, y, z), got {len(coord)}")
        return v

    @validator("fragments")
    @classmethod
    def validate_fragments(cls, v, values):
        """Ensure fragments cover all atoms and don't overlap."""
        if "symbols" in values:
            natoms = len(values["symbols"])
            all_atoms = set()
            for frag in v:
                for atom_idx in frag:
                    if atom_idx < 0 or atom_idx >= natoms:
                        raise ValueError(f"Atom index {atom_idx} out of range [0, {natoms-1}]")
                    if atom_idx in all_atoms:
                        raise ValueError(f"Atom {atom_idx} appears in multiple fragments")
                    all_atoms.add(atom_idx)
            if len(all_atoms) != natoms:
                missing = set(range(natoms)) - all_atoms
                raise ValueError(f"Not all atoms assigned to fragments. Missing: {sorted(missing)}")
        return v


class FileMoleculeSchema(BaseModel):
    """Molecule from external file."""

    file: str = Field(
        ...,
        description="Path to molecule file (XYZ, PDB, or QCSchema JSON).",
    )
    fragments: Optional[List[List[int]]] = Field(
        None,
        description="Fragment definitions. If not specified, will try to auto-detect from file or treat as single fragment.",
    )


class MoleculeSchema(BaseModel):
    """Molecule specification."""

    source: MoleculeSourceEnum = Field(
        MoleculeSourceEnum.inline,
        description="How the molecule is specified.",
    )
    # For inline specification
    inline: Optional[InlineMoleculeSchema] = Field(
        None,
        description="Inline molecule specification. Required if source='inline'.",
    )
    # For file-based specification
    file: Optional[str] = Field(
        None,
        description="Path to molecule file. Required if source='xyz', 'pdb', or 'qcschema'.",
    )
    fragments: Optional[List[List[int]]] = Field(
        None,
        description="Fragment definitions for file-based molecules.",
    )

    @validator("inline")
    @classmethod
    def validate_inline(cls, v, values):
        """Ensure inline is provided when source='inline'."""
        if values.get("source") == MoleculeSourceEnum.inline and v is None:
            raise ValueError("'inline' must be specified when source='inline'")
        return v

    @validator("file")
    @classmethod
    def validate_file(cls, v, values):
        """Ensure file is provided for file-based sources."""
        source = values.get("source")
        if source in [MoleculeSourceEnum.xyz, MoleculeSourceEnum.pdb, MoleculeSourceEnum.qcschema]:
            if v is None:
                raise ValueError(f"'file' must be specified when source='{source.value}'")
        return v


# ==================== Calculation Specification ====================


class SingleLevelCalculationSchema(BaseModel):
    """Single-level calculation specification."""

    driver: DriverEnum = Field(
        DriverEnum.energy,
        description="Type of calculation: energy, gradient, or hessian.",
    )
    method: str = Field(
        ...,
        description="QC method (e.g., 'scf', 'mp2', 'ccsd(t)').",
    )
    basis: str = Field(
        ...,
        description="Basis set (e.g., 'cc-pvdz', 'aug-cc-pvtz').",
    )
    program: str = Field(
        ...,
        description="QC program to use (e.g., 'psi4', 'nwchem', 'cfour').",
    )


class MultiLevelCalculationSchema(BaseModel):
    """Multi-level calculation specification with different methods at different n-body levels."""

    driver: DriverEnum = Field(
        DriverEnum.energy,
        description="Type of calculation: energy, gradient, or hessian.",
    )
    levels: Dict[Union[int, Literal["supersystem"]], Dict[str, str]] = Field(
        ...,
        description="Method specification for each n-body level. Keys are integers (1, 2, 3, ...) or 'supersystem'. "
        "Values are dicts with 'method', 'basis', and 'program'.",
    )

    @validator("levels")
    @classmethod
    def validate_levels(cls, v):
        """Ensure level specifications are valid."""
        for level, spec in v.items():
            # Validate level key
            if level != "supersystem":
                try:
                    level_int = int(level)
                    if level_int < 1:
                        raise ValueError(f"Level must be positive integer, got {level}")
                except (ValueError, TypeError):
                    raise ValueError(f"Level must be integer or 'supersystem', got {level}")

            # Validate spec has required keys
            required_keys = {"method", "basis", "program"}
            if not isinstance(spec, dict):
                raise ValueError(f"Level specification must be dict, got {type(spec)}")
            missing = required_keys - set(spec.keys())
            if missing:
                raise ValueError(f"Level {level} missing required keys: {missing}")

        return v


class CalculationSchema(BaseModel):
    """Calculation specification - either single-level or multi-level."""

    # Use discriminated union pattern
    type: Literal["single", "multi"] = Field(
        "single",
        description="Whether this is a single-level or multi-level calculation.",
    )
    single: Optional[SingleLevelCalculationSchema] = Field(
        None,
        description="Single-level calculation specification.",
    )
    multi: Optional[MultiLevelCalculationSchema] = Field(
        None,
        description="Multi-level calculation specification.",
    )

    @validator("single")
    @classmethod
    def validate_single(cls, v, values):
        """Ensure single is provided when type='single'."""
        if values.get("type") == "single" and v is None:
            raise ValueError("'single' must be specified when type='single'")
        return v

    @validator("multi")
    @classmethod
    def validate_multi(cls, v, values):
        """Ensure multi is provided when type='multi'."""
        if values.get("type") == "multi" and v is None:
            raise ValueError("'multi' must be specified when type='multi'")
        return v


# ==================== Many-Body Options ====================


class BsseSchema(BaseModel):
    """BSSE correction options."""

    type: List[BsseTypeEnum] = Field(
        [BsseTypeEnum.cp],
        description="BSSE correction type(s) to apply. First in list determines returned result.",
    )

    @validator("type", pre=True)
    @classmethod
    def ensure_list(cls, v):
        """Ensure type is a list."""
        if not isinstance(v, list):
            return [v]
        return v


class ManyBodySchema(BaseModel):
    """Many-body expansion options."""

    max_nbody: Optional[int] = Field(
        None,
        description="Maximum n-body level to compute. Default: number of fragments.",
    )
    return_total_data: Optional[bool] = Field(
        None,
        description="Return total energy/gradient/Hessian (True) or interaction data (False). "
        "Default: False for energy, True for gradient/Hessian.",
    )
    supersystem_ie_only: Optional[bool] = Field(
        False,
        description="Only compute supersystem and 1-body for interaction energy. "
        "Skips intermediate n-body calculations. Requires max_nbody=nfragments.",
    )
    embedding_charges: Optional[Dict[int, List[float]]] = Field(
        None,
        description="Point charges for fragments. Keys are 1-based fragment indices. "
        "Requires QCMANYBODY_EMBEDDING_CHARGES environment variable.",
    )


# ==================== Program-Specific Keywords ====================


class ProgramKeywordsSchema(BaseModel):
    """Program-specific keywords."""

    keywords: Optional[Dict[str, Any]] = Field(
        {},
        description="Program-specific keywords (e.g., {'scf__e_convergence': 1e-8}).",
    )


# ==================== Output Options ====================


class OutputSchema(BaseModel):
    """Output configuration."""

    format: OutputFormatEnum = Field(
        OutputFormatEnum.json,
        description="Output format.",
    )
    file: Optional[str] = Field(
        None,
        description="Output file path. If not specified, writes to stdout.",
    )
    pretty: Optional[bool] = Field(
        True,
        description="Pretty-print output (for JSON/YAML).",
    )
    include_timings: Optional[bool] = Field(
        True,
        description="Include timing information in output.",
    )


# ==================== Top-Level Input Schema ====================


class QCManyBodyInput(BaseModel):
    """
    QCManyBody CLI Input File Schema

    This defines the structure of JSON/YAML input files for the QCManyBody CLI.
    """

    schema_name: Literal["qcmanybody_cli_input"] = "qcmanybody_cli_input"
    schema_version: Literal[1] = 1

    molecule: MoleculeSchema = Field(
        ...,
        description="Molecule specification.",
    )
    calculation: CalculationSchema = Field(
        ...,
        description="Calculation specification (single-level or multi-level).",
    )
    bsse: Optional[BsseSchema] = Field(
        BsseSchema(),
        description="BSSE correction options.",
    )
    manybody: Optional[ManyBodySchema] = Field(
        ManyBodySchema(),
        description="Many-body expansion options.",
    )
    program: Optional[ProgramKeywordsSchema] = Field(
        ProgramKeywordsSchema(),
        description="Program-specific keywords.",
    )
    output: Optional[OutputSchema] = Field(
        OutputSchema(),
        description="Output configuration.",
    )

    class Config:
        """Pydantic configuration."""

        # Allow validation of enum string values
        use_enum_values = False
        # Additional config
        extra = "forbid"  # Reject unknown fields
        validate_assignment = True

    def get_driver(self) -> DriverEnum:
        """Get the driver from calculation specification."""
        if self.calculation.type == "single":
            return self.calculation.single.driver
        else:
            return self.calculation.multi.driver

    def get_program(self) -> Optional[str]:
        """Get the program (for single-level calculations only)."""
        if self.calculation.type == "single":
            return self.calculation.single.program
        return None

    def is_multilevel(self) -> bool:
        """Check if this is a multi-level calculation."""
        return self.calculation.type == "multi"
