"""
QCManyBody CLI Run Command

Executes a many-body expansion calculation from an input file.
"""

import json
import logging
from argparse import Namespace
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def handle_run(args: Namespace) -> int:
    """
    Handle the 'run' command.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments

    Returns
    -------
    int
        Exit code (0 for success, non-zero for failure)
    """
    logger.info("Starting QCManyBody calculation")
    logger.info("Input file: %s", args.input)

    # Step 1: Load and validate input file
    try:
        from qcmanybody.cli.input_parser import parse_input_file

        logger.debug("Parsing input file...")
        cli_input = parse_input_file(args.input)
        logger.info(
            "✓ Input file parsed successfully (schema version %s)",
            cli_input.schema_version,
        )

    except Exception as exc:
        logger.error("Failed to parse input file: %s", exc)
        return 1

    # Step 2: Convert to ManyBodyInput
    try:
        from qcmanybody.cli.converter import convert_to_manybody_input

        logger.debug("Converting to ManyBodyInput...")
        mb_input = convert_to_manybody_input(
            cli_input,
            input_file_path=args.input,
        )
        logger.info("✓ Input converted successfully")
        logger.info(
            "  Molecule: %s atoms, %s fragments",
            len(mb_input.molecule.symbols),
            len(mb_input.molecule.fragments),
        )
        logger.info("  Driver: %s", mb_input.specification.driver)
        logger.info(
            "  BSSE types: %s",
            [bt.value for bt in mb_input.specification.keywords.bsse_type],
        )

    except Exception as exc:
        logger.error("Failed to convert input: %s", exc)
        return 1

    # Step 3: Execute calculation via ManyBodyComputer (with parallel support)
    try:
        # Determine execution configuration
        execution_config = cli_input.execution or {}
        parallel = getattr(execution_config, 'parallel', False)
        n_workers = getattr(execution_config, 'n_workers', None)
        executor_type = getattr(execution_config, 'executor_type', 'multiprocessing')
        timeout = getattr(execution_config, 'timeout_per_task', 3600.0)
        max_retries = getattr(execution_config, 'max_retries', 2)

        # Override with command-line arguments if provided
        if hasattr(args, 'parallel') and args.parallel is not None:
            parallel = args.parallel
        if hasattr(args, 'n_workers') and args.n_workers is not None:
            n_workers = args.n_workers

        # Import appropriate computer class
        if parallel:
            from qcmanybody.parallel import ParallelManyBodyComputer, ExecutorConfig
            from qcmanybody.parallel.executors import MultiprocessingExecutor, SequentialExecutor

            logger.info("Initializing ParallelManyBodyComputer...")
            logger.info(f"  Parallel execution: enabled")
            logger.info(f"  Workers: {n_workers if n_workers else 'auto-detect'}")
            logger.info(f"  Executor type: {executor_type}")

            # Create executor config
            config = ExecutorConfig(
                n_workers=n_workers,
                timeout_per_task=timeout,
                max_retries=max_retries,
            )

            # Create executor
            if executor_type == 'sequential':
                executor = SequentialExecutor(config)
            else:  # multiprocessing (default)
                executor = MultiprocessingExecutor(config)

            # Execute with parallel computer
            result = ParallelManyBodyComputer.from_manybodyinput(
                mb_input,
                executor=executor
            )
        else:
            from qcmanybody import ManyBodyComputer

            logger.info("Initializing ManyBodyComputer (sequential mode)...")
            # Note: from_manybodyinput with build_tasks=True (default) executes
            # the calculation and returns a ManyBodyResult, not a ManyBodyComputer
            result = ManyBodyComputer.from_manybodyinput(mb_input)

        logger.info("✓ Calculation completed successfully")

    except ImportError as exc:
        logger.error(
            "Failed to import ManyBodyComputer. Is qcengine installed?"
        )
        logger.error("  Install with: pip install qcengine")
        logger.error("  Error: %s", exc)
        return 1
    except Exception as exc:
        logger.error("Calculation failed: %s", exc)
        logger.debug("Full error:", exc_info=True)
        return 1

    # Determine output preferences (CLI overrides schema)
    output_spec = getattr(cli_input, "output", None)

    if args.format:
        output_format = args.format
    elif output_spec and getattr(output_spec, "format", None):
        output_format = getattr(
            output_spec.format,
            "value",
            output_spec.format,
        )
    else:
        output_format = "json"

    if args.output:
        output_path = Path(args.output)
    elif output_spec and getattr(output_spec, "file", None):
        output_path = Path(output_spec.file)
    else:
        output_path = None

    # Step 4: Format output
    try:
        logger.debug("Formatting output as %s...", output_format)
        output_data = format_output(result, cli_input, output_format)

    except Exception as e:
        logger.error("Failed to format output: %s", e)
        return 1

    # Step 5: Write output
    try:
        if output_path:
            logger.info("Writing output to %s...", output_path)
            output_path.write_text(output_data)
            logger.info("✓ Results written to %s", output_path)
        else:
            # Write to stdout
            print(output_data)

    except Exception as e:
        logger.error("Failed to write output: %s", e)
        return 1

    logger.info("✓ Run command completed successfully")
    return 0


def format_output(result: Any, cli_input: Any, output_format: str) -> str:
    """
    Format the ManyBodyResult for output.

    Parameters
    ----------
    result : ManyBodyResult
        The calculation result
    cli_input : QCManyBodyInput
        The original CLI input
    output_format : str
        Output format: 'json', 'yaml', or 'text'

    Returns
    -------
    str
        Formatted output string
    """
    if output_format == "json":
        return format_json(result, cli_input)
    elif output_format == "yaml":
        return format_yaml(result, cli_input)
    elif output_format == "text":
        return format_text(result, cli_input)
    else:
        raise ValueError(f"Unknown output format: {output_format}")


def format_json(result: Any, cli_input: Any) -> str:
    """Format output as JSON."""
    # Convert result to dictionary
    output = {
        "schema_version": "1.0",
        "input": cli_input.dict(),
        "results": result.dict() if hasattr(result, "dict") else str(result),
    }
    return json.dumps(output, indent=2)


def format_yaml(result: Any, cli_input: Any) -> str:
    """Format output as YAML."""
    try:
        import yaml
    except ImportError:
        logger.warning("PyYAML not installed, falling back to JSON format")
        return format_json(result, cli_input)

    output = {
        "schema_version": "1.0",
        "input": cli_input.dict(),
        "results": result.dict() if hasattr(result, "dict") else str(result),
    }
    return yaml.dump(output, default_flow_style=False, sort_keys=False)


def format_text(result: Any, cli_input: Any) -> str:
    """Format output as human-readable text summary."""
    lines = []
    lines.append("=" * 70)
    lines.append("QCManyBody Calculation Results")
    lines.append("=" * 70)
    lines.append("")

    # Input summary
    lines.append("Input Summary:")
    lines.append("-" * 70)
    mol = cli_input.molecule
    if hasattr(mol, "symbols"):
        lines.append(
            "  Molecule: %s atoms, %s fragments"
            % (len(mol.symbols), len(mol.fragments))
        )
    calc = cli_input.calculation
    driver_summary = "unknown"

    if calc.type == "multi" and calc.multi:
        lines.append("  Calculation: Multi-level")
        for level, spec in calc.multi.levels.items():
            # spec comes from validated schema so contains method/basis/program
            method = spec.get("method", "?")
            basis = spec.get("basis", "?")
            program = spec.get("program", "?")
            lines.append(f"    Level {level}: {method}/{basis} ({program})")
        driver_summary = calc.multi.driver.value
    elif calc.type == "single" and calc.single:
        single = calc.single
        lines.append(
            f"  Method: {single.method}/{single.basis} ({single.program})"
        )
        driver_summary = single.driver.value
    else:
        lines.append("  Calculation details unavailable")

    lines.append(f"  Driver: {driver_summary}")

    bsse_types = []
    if cli_input.bsse and getattr(cli_input.bsse, "type", None):
        bsse_types = [
            bt.value if hasattr(bt, "value") else str(bt)
            for bt in cli_input.bsse.type
        ]

    lines.append(
        "  BSSE: %s"
        % (", ".join(bsse_types) if bsse_types else "default")
    )
    lines.append("")

    # Results summary
    lines.append("Results:")
    lines.append("-" * 70)
    lines.append(f"{result}")
    lines.append("")

    lines.append("=" * 70)
    return "\n".join(lines)
