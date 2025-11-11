"""
QCManyBody CLI Convert Command

Converts input files between JSON and YAML formats.
"""

import json
import logging
from argparse import Namespace
from pathlib import Path

logger = logging.getLogger(__name__)


def handle_convert(args: Namespace) -> int:
    """
    Handle the 'convert' command.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments

    Returns
    -------
    int
        Exit code (0 for success, 1 for failure)
    """
    logger.info("Converting input file")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Output file: {args.output}")

    # Step 1: Parse input file
    try:
        from qcmanybody.cli.input_parser import parse_input_file

        logger.debug("Parsing input file...")
        cli_input = parse_input_file(args.input)
        logger.info(f"✓ Input file parsed successfully")
        logger.info(f"  Format: {'JSON' if args.input.endswith('.json') else 'YAML'}")
        logger.info(f"  Schema version: {cli_input.schema_version}")

    except Exception as e:
        logger.error(f"Failed to parse input file: {e}")
        return 1

    # Step 2: Validate by converting to ManyBodyInput
    try:
        from qcmanybody.cli.converter import convert_to_manybody_input

        logger.debug("Validating input via conversion...")
        mb_input = convert_to_manybody_input(cli_input, input_file_path=args.input)
        logger.info(f"✓ Input is valid")

    except Exception as e:
        logger.error(f"Input validation failed: {e}")
        logger.error(f"  Cannot convert invalid input file")
        return 1

    # Step 3: Determine output format
    output_path = Path(args.output)
    output_ext = output_path.suffix.lower()

    if output_ext == ".json":
        output_format = "json"
    elif output_ext in [".yaml", ".yml"]:
        output_format = "yaml"
    else:
        logger.error(f"Unknown output format: {output_ext}")
        logger.error(f"  Supported formats: .json, .yaml, .yml")
        return 1

    logger.info(f"  Output format: {output_format.upper()}")

    # Step 4: Convert to output format
    try:
        if output_format == "json":
            output_data = format_as_json(cli_input)
        else:  # yaml
            output_data = format_as_yaml(cli_input)

    except Exception as e:
        logger.error(f"Failed to format output: {e}")
        return 1

    # Step 5: Write output file
    try:
        logger.debug(f"Writing to {output_path}...")
        output_path.write_text(output_data)
        logger.info(f"✓ Successfully converted to {output_path}")

        # Show conversion summary
        input_size = Path(args.input).stat().st_size
        output_size = output_path.stat().st_size
        print()
        print(f"Conversion Summary:")
        print(f"  Input:  {args.input} ({input_size} bytes)")
        print(f"  Output: {args.output} ({output_size} bytes)")
        print(f"  ✓ Conversion successful")
        print()

    except Exception as e:
        logger.error(f"Failed to write output file: {e}")
        return 1

    return 0


def format_as_json(cli_input) -> str:
    """
    Format CLI input as JSON.

    Parameters
    ----------
    cli_input : QCManyBodyInput
        Parsed CLI input

    Returns
    -------
    str
        JSON formatted string
    """
    # Use Pydantic's json() method which properly serializes enums
    json_str = cli_input.json(exclude_none=True, indent=2)
    return json_str + "\n"


def format_as_yaml(cli_input) -> str:
    """
    Format CLI input as YAML.

    Parameters
    ----------
    cli_input : QCManyBodyInput
        Parsed CLI input

    Returns
    -------
    str
        YAML formatted string
    """
    try:
        import yaml
    except ImportError:
        raise RuntimeError(
            "PyYAML is required for YAML output but is not installed.\n" "Install with: pip install pyyaml"
        )

    # Convert via JSON to ensure proper enum serialization, then to Python dict
    json_str = cli_input.json(exclude_none=True)
    data = json.loads(json_str)
    return yaml.dump(data, default_flow_style=False, sort_keys=False, allow_unicode=True)
