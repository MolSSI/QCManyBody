"""
QCManyBody CLI Input Converter

Converts CLI input schema to internal ManyBodyInput format.
"""

import logging
from typing import Any, Dict, Optional

from qcelemental.models import Molecule

from qcmanybody.cli.molecule_loader import MoleculeLoadError, load_molecule
from qcmanybody.cli.schemas.input_schema import QCManyBodyInput
from qcmanybody.models import (
    AtomicSpecification,
    BsseEnum,
    ManyBodyInput,
    ManyBodyKeywords,
    ManyBodyProtocols,
    ManyBodySpecification,
)

logger = logging.getLogger(__name__)


class ConversionError(Exception):
    """Exception raised when conversion fails."""

    pass


def convert_bsse_types(bsse_types: list) -> list:
    """
    Convert CLI BSSE types to internal BsseEnum.

    Parameters
    ----------
    bsse_types : list
        List of BSSE type strings from CLI input

    Returns
    -------
    list
        List of BsseEnum values
    """
    result = []
    for bt in bsse_types:
        # Convert to BsseEnum
        try:
            result.append(BsseEnum[bt.value if hasattr(bt, "value") else bt])
        except KeyError:
            # Already should be validated, but just in case
            raise ConversionError(f"Invalid BSSE type: {bt}")

    return result


def create_atomic_specification(
    method: str, basis: str, program: str, driver: str, keywords: Optional[Dict[str, Any]] = None
) -> AtomicSpecification:
    """
    Create AtomicSpecification from CLI parameters.

    Parameters
    ----------
    method : str
        QC method
    basis : str
        Basis set
    program : str
        QC program
    driver : str
        Driver type
    keywords : Optional[Dict[str, Any]]
        Program-specific keywords

    Returns
    -------
    AtomicSpecification
        Atomic specification object
    """
    if keywords is None:
        keywords = {}

    return AtomicSpecification(
        program=program,
        driver=driver,
        model={"method": method, "basis": basis},
        keywords=keywords,
        protocols={"stdout": False},  # Don't store stdout by default for efficiency
    )


def convert_single_level(cli_input: QCManyBodyInput) -> Dict[str, AtomicSpecification]:
    """
    Convert single-level calculation specification.

    Parameters
    ----------
    cli_input : QCManyBodyInput
        CLI input object

    Returns
    -------
    Dict[str, AtomicSpecification]
        Specification dictionary for ManyBodyInput
    """
    single = cli_input.calculation.single

    # Get program keywords
    keywords = cli_input.program.keywords if cli_input.program else {}

    # Create atomic specification
    atomic_spec = create_atomic_specification(
        method=single.method,
        basis=single.basis,
        program=single.program,
        driver=single.driver.value,
        keywords=keywords,
    )

    # Return as dict with auto key
    return {"(auto)": atomic_spec}


def convert_multi_level(cli_input: QCManyBodyInput) -> Dict[str, AtomicSpecification]:
    """
    Convert multi-level calculation specification.

    Parameters
    ----------
    cli_input : QCManyBodyInput
        CLI input object

    Returns
    -------
    Dict[str, AtomicSpecification]
        Specification dictionary for ManyBodyInput with level keys

    Notes
    -----
    Multi-level format in CLI:
    {
      "levels": {
        "1": {"method": "ccsd(t)", "basis": "cc-pvtz", "program": "nwchem"},
        "2": {"method": "mp2", "basis": "cc-pvdz", "program": "psi4"},
        "supersystem": {"method": "scf", "basis": "cc-pvdz", "program": "psi4"}
      }
    }

    Internal format uses method/basis strings as keys for the levels field in keywords,
    and program-specific specs in the specification field.
    """
    multi = cli_input.calculation.multi
    driver = multi.driver.value

    # Get base keywords (same for all levels)
    base_keywords = cli_input.program.keywords if cli_input.program else {}

    # Build levels dict for ManyBodyKeywords
    # This maps level number -> "method/basis" string
    levels_dict = {}
    specifications_dict = {}

    for level_key, level_spec in multi.levels.items():
        method = level_spec["method"]
        basis = level_spec["basis"]
        program = level_spec["program"]

        # The internal format uses "method/basis" as the specification key
        method_basis_key = f"{method}/{basis}"

        # Store in levels dict
        if level_key == "supersystem":
            levels_dict["supersystem"] = method_basis_key
        else:
            levels_dict[int(level_key)] = method_basis_key

        # Create atomic specification for this method/basis/program combination
        # Only create if we haven't seen this combination before
        if method_basis_key not in specifications_dict:
            # Get program-specific keywords for this level
            # For now, use the same keywords for all levels
            # Could be enhanced to support per-level keywords
            atomic_spec = create_atomic_specification(
                method=method,
                basis=basis,
                program=program,
                driver=driver,
                keywords=base_keywords,
            )
            specifications_dict[method_basis_key] = atomic_spec

    # For multi-level, we need to return the specifications dict
    # The levels dict goes into ManyBodyKeywords
    return specifications_dict, levels_dict


def create_manybody_keywords(cli_input: QCManyBodyInput, levels_dict: Optional[Dict] = None) -> ManyBodyKeywords:
    """
    Create ManyBodyKeywords from CLI input.

    Parameters
    ----------
    cli_input : QCManyBodyInput
        CLI input object
    levels_dict : Optional[Dict]
        Levels dictionary for multi-level calculations

    Returns
    -------
    ManyBodyKeywords
        ManyBody keywords object
    """
    # Convert BSSE types
    bsse_types = convert_bsse_types(cli_input.bsse.type)

    # Build keywords dict
    kwargs = {
        "bsse_type": bsse_types,
    }

    # Add optional fields from manybody section
    if cli_input.manybody:
        if cli_input.manybody.max_nbody is not None:
            kwargs["max_nbody"] = cli_input.manybody.max_nbody

        if cli_input.manybody.return_total_data is not None:
            kwargs["return_total_data"] = cli_input.manybody.return_total_data

        if cli_input.manybody.supersystem_ie_only is not None:
            kwargs["supersystem_ie_only"] = cli_input.manybody.supersystem_ie_only

        if cli_input.manybody.embedding_charges is not None:
            kwargs["embedding_charges"] = cli_input.manybody.embedding_charges

    # Add levels for multi-level calculations
    if levels_dict is not None:
        kwargs["levels"] = levels_dict

    return ManyBodyKeywords(**kwargs)


def create_manybody_protocols(cli_input: QCManyBodyInput) -> ManyBodyProtocols:
    """
    Create ManyBodyProtocols from CLI input.

    Parameters
    ----------
    cli_input : QCManyBodyInput
        CLI input object

    Returns
    -------
    ManyBodyProtocols
        ManyBody protocols object
    """
    # For now, use default: don't store component results (for efficiency)
    # Could be made configurable via CLI options
    return ManyBodyProtocols(component_results="none")


def convert_to_manybody_input(cli_input: QCManyBodyInput, input_file_path: Optional[str] = None) -> ManyBodyInput:
    """
    Convert CLI input schema to internal ManyBodyInput.

    This is the main conversion function.

    Parameters
    ----------
    cli_input : QCManyBodyInput
        Validated CLI input object
    input_file_path : Optional[str]
        Path to the input file (for resolving relative paths in molecule files).
        If provided, relative paths will be resolved relative to this file's directory.

    Returns
    -------
    ManyBodyInput
        Internal ManyBodyInput object ready for computation

    Raises
    ------
    ConversionError
        If conversion fails
    MoleculeLoadError
        If molecule loading fails

    Examples
    --------
    >>> from qcmanybody.cli.input_parser import parse_input_file
    >>> from qcmanybody.cli.converter import convert_to_manybody_input
    >>> cli_input = parse_input_file("my_calc.json")
    >>> mbe_input = convert_to_manybody_input(cli_input, "my_calc.json")
    >>> # Now use mbe_input with ManyBodyComputer
    """
    logger.info("Converting CLI input to ManyBodyInput")

    # Resolve molecule file path if necessary
    mol_schema = cli_input.molecule
    if input_file_path and mol_schema.file:
        from pathlib import Path

        input_dir = Path(input_file_path).parent
        mol_file = Path(mol_schema.file)

        # If molecule file path is relative, make it relative to input file
        if not mol_file.is_absolute():
            resolved_path = input_dir / mol_file
            logger.debug(f"Resolving relative path: {mol_schema.file} -> {resolved_path}")
            # Create a modified schema with resolved path
            from qcmanybody.cli.schemas.input_schema import MoleculeSchema

            mol_schema = MoleculeSchema(**{**mol_schema.dict(), "file": str(resolved_path)})

    # Load molecule
    try:
        molecule = load_molecule(mol_schema)
    except MoleculeLoadError as e:
        raise ConversionError(f"Failed to load molecule: {e}") from e

    # Get driver
    driver = cli_input.get_driver().value

    # Convert calculation specification
    if cli_input.calculation.type == "single":
        specifications = convert_single_level(cli_input)
        levels_dict = None
    else:
        specifications, levels_dict = convert_multi_level(cli_input)

    # Create ManyBody keywords
    mb_keywords = create_manybody_keywords(cli_input, levels_dict)

    # Create ManyBody protocols
    mb_protocols = create_manybody_protocols(cli_input)

    # Create ManyBodySpecification
    mb_spec = ManyBodySpecification(
        keywords=mb_keywords,
        driver=driver,
        specification=specifications,
        protocols=mb_protocols,
    )

    # Create ManyBodyInput
    mb_input = ManyBodyInput(
        molecule=molecule,
        specification=mb_spec,
    )

    logger.info("Successfully converted to ManyBodyInput")
    logger.debug(f"  Molecule: {len(molecule.symbols)} atoms, {len(molecule.fragments)} fragments")
    logger.debug(f"  Driver: {driver}")
    logger.debug(f"  BSSE types: {[bt.value for bt in mb_keywords.bsse_type]}")
    if mb_keywords.max_nbody:
        logger.debug(f"  Max n-body: {mb_keywords.max_nbody}")

    return mb_input
