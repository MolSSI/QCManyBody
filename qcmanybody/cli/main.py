"""
QCManyBody CLI Main Entry Point

This module provides the main entry point for the QCManyBody command-line interface.
"""

import argparse
import logging
import sys
from typing import Optional

# Try to get version from package
try:
    from qcmanybody import __version__
except ImportError:
    __version__ = "unknown"


def setup_logging(verbosity: int, log_file: Optional[str] = None) -> None:
    """
    Configure logging based on verbosity level.

    Parameters
    ----------
    verbosity : int
        Verbosity level (0=WARNING, 1=INFO, 2+=DEBUG)
    log_file : str, optional
        Path to log file. If None, log to stderr only.
    """
    levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    level = levels[min(verbosity, len(levels) - 1)]

    # Configure logging format
    log_format = "%(levelname)s: %(message)s"
    if verbosity >= 2:
        log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Set up handlers
    handlers = []

    # Console handler (stderr)
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(logging.Formatter(log_format))
    handlers.append(console_handler)

    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)  # Always log everything to file
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        handlers.append(file_handler)

    # Configure root logger
    logging.basicConfig(level=level, format=log_format, handlers=handlers, force=True)


def create_parser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser for the CLI.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser
    """
    # Create main parser
    parser = argparse.ArgumentParser(
        prog="qcmanybody",
        description="QCManyBody: Many-body expansion calculations for quantum chemistry",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  qcmanybody run input.json -o results.json
  qcmanybody plan input.json --show-tasks
  qcmanybody validate input.json
  qcmanybody run input.json --format yaml -o results.yaml

For more information, visit: https://github.com/MolSSI/QCManyBody
        """,
    )

    # Global options
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (can be repeated: -v, -vv, -vvv)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress non-essential output")

    # Create subparsers for commands
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # ===== RUN COMMAND =====
    run_parser = subparsers.add_parser(
        "run",
        help="Run a QCManyBody calculation",
        description="Execute a many-body expansion calculation from an input file",
    )
    run_parser.add_argument("input", help="Input file path (JSON or YAML)")
    run_parser.add_argument("-o", "--output", help="Output file path (defaults to schema setting or stdout)")
    run_parser.add_argument(
        "--format",
        choices=["json", "yaml", "text"],
        default=None,
        help="Override output format (json, yaml, text)",
    )

    # Parallel execution options
    parallel_group = run_parser.add_argument_group("parallel execution options")
    parallel_group.add_argument(
        "--parallel",
        action="store_true",
        default=None,
        help="Enable parallel execution (overrides input file setting)",
    )
    parallel_group.add_argument(
        "--no-parallel",
        action="store_false",
        dest="parallel",
        help="Disable parallel execution (force sequential mode)",
    )
    parallel_group.add_argument(
        "--n-workers",
        type=int,
        default=None,
        metavar="N",
        help="Number of parallel workers (overrides input file setting, auto-detect if not specified)",
    )
    parallel_group.add_argument(
        "--schedule-strategy",
        choices=["priority_first", "cost_first", "nbody_first", "dependency_aware", "fifo"],
        default=None,
        help="Task scheduling strategy for parallel execution (default: priority_first)",
    )
    parallel_group.add_argument(
        "--checkpoint-file",
        type=str,
        default=None,
        metavar="FILE",
        help="Path to checkpoint file for saving/resuming progress",
    )
    parallel_group.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Resume from checkpoint file if it exists (requires --checkpoint-file)",
    )

    # Logging options
    log_group = run_parser.add_argument_group("logging options")
    log_group.add_argument("--log", help="Log file path")

    # ===== PLAN COMMAND =====
    plan_parser = subparsers.add_parser(
        "plan",
        help="Show execution plan without running",
        description="Display the tasks that would be executed for a calculation",
    )
    plan_parser.add_argument("input", help="Input file path (JSON or YAML)")
    plan_parser.add_argument("--show-tasks", action="store_true", help="List all individual tasks")
    plan_parser.add_argument("--log", help="Log file path")

    # ===== VALIDATE COMMAND =====
    validate_parser = subparsers.add_parser(
        "validate",
        help="Validate input file",
        description="Check input file for errors without running calculation",
    )
    validate_parser.add_argument(
        "input", nargs="?", help="Input file path (JSON or YAML, not required if --show-schema is used)"
    )
    validate_parser.add_argument("--strict", action="store_true", help="Enable strict validation mode")
    validate_parser.add_argument("--show-schema", action="store_true", help="Display expected schema")
    validate_parser.add_argument("--log", help="Log file path")

    # ===== CONVERT COMMAND =====
    convert_parser = subparsers.add_parser(
        "convert", help="Convert between input formats", description="Convert input files between JSON and YAML formats"
    )
    convert_parser.add_argument("input", help="Input file path")
    convert_parser.add_argument("output", help="Output file path")
    convert_parser.add_argument("--log", help="Log file path")

    return parser


def main(argv: Optional[list] = None) -> int:
    """
    Main entry point for the QCManyBody CLI.

    Parameters
    ----------
    argv : list, optional
        Command-line arguments. If None, uses sys.argv[1:]

    Returns
    -------
    int
        Exit code (0 for success, non-zero for failure)
    """
    # Create parser
    parser = create_parser()

    # Parse arguments
    args = parser.parse_args(argv)

    # Set up logging
    verbosity = 0 if args.quiet else args.verbose
    log_file = getattr(args, "log", None)
    setup_logging(verbosity, log_file)

    logger = logging.getLogger(__name__)
    logger.debug(f"QCManyBody CLI v{__version__}")
    logger.debug(f"Command: {args.command}")
    logger.debug(f"Arguments: {args}")

    # Route to appropriate command handler
    try:
        if args.command == "run":
            from qcmanybody.cli.commands.run import handle_run

            return handle_run(args)

        elif args.command == "plan":
            from qcmanybody.cli.commands.plan import handle_plan

            return handle_plan(args)

        elif args.command == "validate":
            from qcmanybody.cli.commands.validate import handle_validate

            return handle_validate(args)

        elif args.command == "convert":
            from qcmanybody.cli.commands.convert import handle_convert

            return handle_convert(args)

        else:
            logger.error(f"Unknown command: {args.command}")
            parser.print_help()
            return 1

    except KeyboardInterrupt:
        logger.warning("Interrupted by user")
        return 130  # Standard exit code for SIGINT

    except Exception as e:
        logger.error(f"Error: {e}", exc_info=verbosity >= 2)
        return 1


if __name__ == "__main__":
    sys.exit(main())
