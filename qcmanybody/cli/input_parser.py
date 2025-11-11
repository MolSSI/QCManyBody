"""
QCManyBody CLI Input File Parser

Handles loading and validating input files (JSON/YAML).
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

from pydantic.v1 import ValidationError

from qcmanybody.cli.schemas.input_schema import QCManyBodyInput

logger = logging.getLogger(__name__)

# Check for YAML support
try:
    import yaml

    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False
    logger.debug("PyYAML not available - YAML input files not supported")


class InputParseError(Exception):
    """Exception raised when input file parsing fails."""

    pass


class InputValidationError(Exception):
    """Exception raised when input validation fails."""

    pass


def detect_file_format(file_path: Path) -> str:
    """
    Detect input file format from extension.

    Parameters
    ----------
    file_path : Path
        Path to input file

    Returns
    -------
    str
        File format: 'json' or 'yaml'

    Raises
    ------
    InputParseError
        If file format cannot be determined
    """
    suffix = file_path.suffix.lower()

    if suffix in [".json"]:
        return "json"
    elif suffix in [".yaml", ".yml"]:
        return "yaml"
    else:
        raise InputParseError(
            f"Cannot determine file format from extension '{suffix}'. " f"Supported extensions: .json, .yaml, .yml"
        )


def load_json_file(file_path: Path) -> Dict[str, Any]:
    """
    Load JSON file.

    Parameters
    ----------
    file_path : Path
        Path to JSON file

    Returns
    -------
    Dict[str, Any]
        Parsed JSON data

    Raises
    ------
    InputParseError
        If JSON parsing fails
    """
    try:
        with open(file_path, "r") as f:
            data = json.load(f)
        logger.debug(f"Successfully loaded JSON from {file_path}")
        return data
    except json.JSONDecodeError as e:
        raise InputParseError(
            f"Failed to parse JSON file {file_path}:\n"
            f"  Line {e.lineno}, Column {e.colno}: {e.msg}\n"
            f"  Check for missing commas, brackets, or quotes."
        ) from e
    except IOError as e:
        raise InputParseError(f"Failed to read file {file_path}: {e}") from e


def load_yaml_file(file_path: Path) -> Dict[str, Any]:
    """
    Load YAML file.

    Parameters
    ----------
    file_path : Path
        Path to YAML file

    Returns
    -------
    Dict[str, Any]
        Parsed YAML data

    Raises
    ------
    InputParseError
        If YAML parsing fails or PyYAML not available
    """
    if not YAML_AVAILABLE:
        raise InputParseError(
            f"Cannot load YAML file {file_path}: PyYAML is not installed.\n"
            f"Install with: pip install qcmanybody[cli]\n"
            f"Or use JSON format instead (supported without additional dependencies)."
        )

    try:
        with open(file_path, "r") as f:
            data = yaml.safe_load(f)
        logger.debug(f"Successfully loaded YAML from {file_path}")
        return data
    except yaml.YAMLError as e:
        # Format YAML error message
        error_msg = f"Failed to parse YAML file {file_path}:\n"
        if hasattr(e, "problem_mark"):
            mark = e.problem_mark
            error_msg += f"  Line {mark.line + 1}, Column {mark.column + 1}\n"
        if hasattr(e, "problem"):
            error_msg += f"  {e.problem}\n"
        if hasattr(e, "context"):
            error_msg += f"  {e.context}\n"
        raise InputParseError(error_msg) from e
    except IOError as e:
        raise InputParseError(f"Failed to read file {file_path}: {e}") from e


def load_input_file(file_path: str) -> Dict[str, Any]:
    """
    Load input file (auto-detects JSON or YAML).

    Parameters
    ----------
    file_path : str
        Path to input file

    Returns
    -------
    Dict[str, Any]
        Parsed input data

    Raises
    ------
    InputParseError
        If file cannot be loaded or parsed
    """
    path = Path(file_path)

    if not path.exists():
        raise InputParseError(f"Input file not found: {file_path}")

    if not path.is_file():
        raise InputParseError(f"Input path is not a file: {file_path}")

    # Detect format and load
    file_format = detect_file_format(path)
    logger.info(f"Loading {file_format.upper()} input file: {file_path}")

    if file_format == "json":
        return load_json_file(path)
    elif file_format == "yaml":
        return load_yaml_file(path)
    else:
        raise InputParseError(f"Unsupported file format: {file_format}")


def format_validation_error(error: ValidationError, file_path: str) -> str:
    """
    Format Pydantic validation error into user-friendly message.

    Parameters
    ----------
    error : ValidationError
        Pydantic validation error
    file_path : str
        Path to input file (for context)

    Returns
    -------
    str
        Formatted error message
    """
    msg = f"Validation failed for input file: {file_path}\n\n"
    msg += "Errors found:\n"

    for i, err in enumerate(error.errors(), 1):
        # Get field location
        loc = " -> ".join(str(x) for x in err["loc"])

        # Get error message
        err_msg = err["msg"]

        # Get error type
        err_type = err["type"]

        msg += f"\n{i}. Field: {loc}\n"
        msg += f"   Error: {err_msg}\n"

        # Add helpful context for common errors
        if err_type == "value_error.missing":
            msg += f"   Hint: This required field is missing. Add it to your input file.\n"
        elif err_type == "type_error":
            msg += f"   Hint: Check the data type of this field.\n"
        elif err_type == "value_error.const":
            msg += f"   Hint: This field has a fixed value that cannot be changed.\n"

    msg += "\nFor input file format documentation, see: [documentation URL]"

    return msg


def validate_input(data: Dict[str, Any], file_path: str) -> QCManyBodyInput:
    """
    Validate input data against schema.

    Parameters
    ----------
    data : Dict[str, Any]
        Input data to validate
    file_path : str
        Path to input file (for error messages)

    Returns
    -------
    QCManyBodyInput
        Validated input object

    Raises
    ------
    InputValidationError
        If validation fails
    """
    try:
        validated = QCManyBodyInput(**data)
        logger.info("Input validation successful")
        return validated
    except ValidationError as e:
        error_msg = format_validation_error(e, file_path)
        raise InputValidationError(error_msg) from e


def parse_input_file(file_path: str) -> QCManyBodyInput:
    """
    Load and validate input file.

    This is the main entry point for input file parsing.

    Parameters
    ----------
    file_path : str
        Path to input file (JSON or YAML)

    Returns
    -------
    QCManyBodyInput
        Validated input object

    Raises
    ------
    InputParseError
        If file cannot be loaded or parsed
    InputValidationError
        If validation fails

    Examples
    --------
    >>> input_data = parse_input_file("my_calculation.json")
    >>> print(f"Driver: {input_data.get_driver()}")
    >>> print(f"BSSE type: {input_data.bsse.type}")
    """
    # Load file
    try:
        data = load_input_file(file_path)
    except InputParseError:
        raise  # Re-raise with original message

    # Validate against schema
    try:
        validated = validate_input(data, file_path)
    except InputValidationError:
        raise  # Re-raise with original message

    return validated


def get_schema_dict() -> Dict[str, Any]:
    """
    Get the JSON schema for the input format.

    Useful for generating documentation or validation tools.

    Returns
    -------
    Dict[str, Any]
        JSON schema
    """
    return QCManyBodyInput.schema()


def get_schema_json(indent: int = 2) -> str:
    """
    Get the JSON schema as a formatted JSON string.

    Parameters
    ----------
    indent : int
        Indentation level for pretty-printing

    Returns
    -------
    str
        JSON schema as string
    """
    schema = get_schema_dict()
    return json.dumps(schema, indent=indent)
