"""
QCManyBody CLI Validate Command

Validates an input file without executing the calculation.
"""

import json
import logging
from argparse import Namespace

logger = logging.getLogger(__name__)


def handle_validate(args: Namespace) -> int:
    """
    Handle the 'validate' command.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments

    Returns
    -------
    int
        Exit code (0 for success, 1 for validation failure)
    """
    # If --show-schema was requested, display schema and exit
    if args.show_schema:
        return show_schema()

    # Check if input file was provided
    if not args.input:
        print("Error: input file is required unless --show-schema is used")
        print("Usage: qcmanybody validate <input-file>")
        print("       qcmanybody validate --show-schema")
        return 1

    logger.info("Validating input file")
    logger.info(f"Input file: {args.input}")

    validation_errors = []
    validation_warnings = []

    # Step 1: Parse and validate input file
    try:
        from qcmanybody.cli.input_parser import parse_input_file

        logger.debug("Parsing input file...")
        cli_input = parse_input_file(args.input)
        print(f"✓ Input file schema is valid")
        print(f"  Schema version: {cli_input.schema_version}")

    except Exception as e:
        validation_errors.append(f"Input parsing failed: {e}")
        print(f"✗ Input file validation FAILED")
        print(f"  Error: {e}")
        return 1

    # Step 2: Validate molecule
    try:
        print(f"\n✓ Molecule specification is valid")
        mol = cli_input.molecule
        if mol.inline:
            print(f"  Type: Inline molecule")
            print(f"  Atoms: {len(mol.inline.symbols)}")
            print(f"  Fragments: {len(mol.inline.fragments)}")
            # Check that all atoms are covered by fragments
            all_atoms = set()
            for frag in mol.inline.fragments:
                all_atoms.update(frag)
            if len(all_atoms) != len(mol.inline.symbols):
                validation_warnings.append(
                    f"Not all atoms are covered by fragments. "
                    f"Atoms: {len(mol.inline.symbols)}, Fragments cover: {len(all_atoms)}"
                )
        elif mol.source.value == "file" and mol.file:
            print(f"  Type: From file")
            print(f"  File: {mol.file.file}")
            if mol.file.fragments:
                print(f"  Fragments: {len(mol.file.fragments)}")

    except Exception as e:
        validation_errors.append(f"Molecule validation failed: {e}")

    # Step 3: Validate calculation settings
    try:
        print(f"\n✓ Calculation settings are valid")
        calc = cli_input.calculation
        if calc.single:
            print(f"  Type: Single-level")
            print(f"  Driver: {calc.single.driver}")
            print(f"  Method: {calc.single.method}")
            print(f"  Basis: {calc.single.basis}")
            print(f"  Program: {calc.single.program}")
        elif calc.multi:
            print(f"  Type: Multi-level")
            print(f"  Levels: {len(calc.multi.levels)}")
            for level, spec in calc.multi.levels.items():
                print(f"    Level {level}: {spec}")

    except Exception as e:
        validation_errors.append(f"Calculation validation failed: {e}")

    # Step 4: Validate BSSE settings
    try:
        if cli_input.bsse:
            print(f"\n✓ BSSE settings are valid")
            bsse_types = cli_input.bsse.type
            if isinstance(bsse_types, list):
                print(f"  Types: {[t.value if hasattr(t, 'value') else str(t) for t in bsse_types]}")
            else:
                print(f"  Type: {bsse_types.value if hasattr(bsse_types, 'value') else str(bsse_types)}")

    except Exception as e:
        validation_errors.append(f"BSSE validation failed: {e}")

    # Step 5: Validate many-body settings
    try:
        if cli_input.manybody:
            print(f"\n✓ Many-body settings are valid")
            print(f"  Max n-body: {cli_input.manybody.max_nbody}")
            print(f"  Return total data: {cli_input.manybody.return_total_data}")

    except Exception as e:
        validation_errors.append(f"Many-body validation failed: {e}")

    # Step 6: Try converting to ManyBodyInput (comprehensive validation)
    try:
        from qcmanybody.cli.converter import convert_to_manybody_input

        logger.debug("Testing conversion to ManyBodyInput...")
        mb_input = convert_to_manybody_input(cli_input, input_file_path=args.input)
        print(f"\n✓ Successfully converts to ManyBodyInput")
        print(f"  Can be used with ManyBodyComputer")

    except Exception as e:
        validation_errors.append(f"Conversion to ManyBodyInput failed: {e}")
        if args.strict:
            print(f"\n✗ Strict validation FAILED")
            print(f"  Error: {e}")
            return 1

    # Display validation summary
    print("\n" + "=" * 70)
    if validation_errors:
        print("✗ VALIDATION FAILED")
        print(f"\nErrors ({len(validation_errors)}):")
        for i, error in enumerate(validation_errors, 1):
            print(f"  {i}. {error}")
        if validation_warnings:
            print(f"\nWarnings ({len(validation_warnings)}):")
            for i, warning in enumerate(validation_warnings, 1):
                print(f"  {i}. {warning}")
        print("=" * 70)
        return 1
    else:
        print("✓ VALIDATION SUCCESSFUL")
        if validation_warnings:
            print(f"\nWarnings ({len(validation_warnings)}):")
            for i, warning in enumerate(validation_warnings, 1):
                print(f"  {i}. {warning}")
        print("\nThe input file is valid and ready to use.")
        print("=" * 70)
        return 0


def show_schema() -> int:
    """Display the input file schema."""
    from qcmanybody.cli.schemas.input_schema import QCManyBodyInput

    print("QCManyBody CLI Input Schema")
    print("=" * 70)
    print()
    print("The input file should be a JSON or YAML file with the following structure:")
    print()

    # Get schema dict
    schema = QCManyBodyInput.schema()

    # Pretty print schema
    print(json.dumps(schema, indent=2))
    print()
    print("=" * 70)
    print()
    print("For detailed documentation and examples, see:")
    print("  cli_addition/INPUT_FILE_SPEC.md")
    print("  cli_addition/EXAMPLES.md")
    print("  examples/cli/")
    print()

    return 0
