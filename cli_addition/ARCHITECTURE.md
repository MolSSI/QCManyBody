# QCManyBody CLI - Technical Architecture

## Architecture Overview

The CLI will be implemented as a thin layer on top of the existing QCManyBody library, following the principle of minimal invasiveness to the core codebase.

```
┌─────────────────────────────────────────────────────────┐
│                     User Interface                       │
│                   (Terminal/Shell)                       │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│                   CLI Entry Point                        │
│              qcmanybody/cli/main.py                     │
│  - Argument parsing (argparse/click/typer)              │
│  - Command routing                                       │
│  - Global configuration                                  │
└────────────────────┬────────────────────────────────────┘
                     │
        ┌────────────┼────────────┬──────────────┐
        ▼            ▼            ▼              ▼
   ┌─────────┐ ┌─────────┐ ┌──────────┐  ┌───────────┐
   │   run   │ │  plan   │ │ validate │  │  convert  │
   │ command │ │ command │ │ command  │  │  command  │
   └────┬────┘ └────┬────┘ └────┬─────┘  └─────┬─────┘
        │           │            │              │
        └───────────┴────────────┴──────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│              Input File Processing Layer                │
│            qcmanybody/cli/input_parser.py               │
│  - YAML/JSON parsing                                     │
│  - Schema validation                                     │
│  - Molecule loading (XYZ, file refs, etc.)              │
│  - Conversion to ManyBodyInput                          │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│              ManyBodyComputer Integration               │
│         Uses existing ManyBodyComputer API              │
│  - No modifications to existing code                    │
│  - Uses ManyBodyInput as interface                      │
│  - Leverages QCEngine for QC program execution         │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│              Output Processing Layer                     │
│            qcmanybody/cli/output_writer.py              │
│  - Format conversion (JSON, YAML, HDF5, text)           │
│  - Summary generation                                    │
│  - Logging and reporting                                │
└─────────────────────────────────────────────────────────┘
```

## Module Structure

```
qcmanybody/
├── cli/
│   ├── __init__.py              # Package init, exports main CLI
│   ├── main.py                  # Entry point, command routing
│   ├── commands/
│   │   ├── __init__.py
│   │   ├── run.py               # Run calculation command
│   │   ├── plan.py              # Show execution plan
│   │   ├── validate.py          # Validate input file
│   │   └── convert.py           # Convert input formats
│   ├── input_parser.py          # Input file parsing and validation
│   ├── molecule_loader.py       # Molecule loading from various formats
│   ├── converter.py             # CLI input → ManyBodyInput conversion
│   ├── output_writer.py         # Output formatting and writing
│   └── schemas/
│       ├── __init__.py
│       └── input_schema.py      # Input file schema definitions
├── core.py                       # Existing ManyBodyCore
├── computer.py                   # Existing ManyBodyComputer
└── ...

examples/
└── cli/
    ├── 01_basic_energy.json
    ├── 01_basic_energy.yaml
    ├── 02_gradient.json
    ├── 03_multilevel.json
    └── 04_from_xyz.json
```

## Component Details

### 1. CLI Entry Point (`main.py`)

**Responsibilities:**
- Parse command-line arguments
- Route to appropriate command handlers
- Handle global options (verbosity, config files, etc.)
- Set up logging

**Implementation with argparse:**

```python
import argparse
import sys
from qcmanybody import __version__

def setup_logging(verbosity):
    """Configure logging based on verbosity level."""
    import logging
    levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    level = levels[min(verbosity, len(levels)-1)]
    logging.basicConfig(
        level=level,
        format='%(levelname)s: %(message)s'
    )

def main():
    """Main CLI entry point."""
    # Create main parser
    parser = argparse.ArgumentParser(
        prog='qcmanybody',
        description='QCManyBody: Many-body expansion calculations'
    )

    # Global options
    parser.add_argument('--version', action='version',
                       version=f'%(prog)s {__version__}')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                       help='Increase verbosity (can be repeated)')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Suppress non-essential output')

    # Create subparsers for commands
    subparsers = parser.add_subparsers(dest='command', required=True,
                                       help='Available commands')

    # Run command
    run_parser = subparsers.add_parser('run',
                                       help='Run a QCManyBody calculation')
    run_parser.add_argument('input', help='Input file path (JSON or YAML)')
    run_parser.add_argument('-o', '--output', help='Output file path')
    run_parser.add_argument('--format',
                          choices=['json', 'yaml', 'text'],
                          default='json',
                          help='Output format (default: json)')
    run_parser.add_argument('--log', help='Log file path')

    # Plan command
    plan_parser = subparsers.add_parser('plan',
                                        help='Show execution plan')
    plan_parser.add_argument('input', help='Input file path')
    plan_parser.add_argument('--show-tasks', action='store_true',
                           help='List all tasks')

    # Validate command
    validate_parser = subparsers.add_parser('validate',
                                           help='Validate input file')
    validate_parser.add_argument('input', help='Input file path')
    validate_parser.add_argument('--strict', action='store_true',
                               help='Enable strict validation')

    # Convert command
    convert_parser = subparsers.add_parser('convert',
                                          help='Convert between formats')
    convert_parser.add_argument('input', help='Input file path')
    convert_parser.add_argument('output', help='Output file path')

    # Parse arguments
    args = parser.parse_args()

    # Set up logging
    verbosity = 0 if args.quiet else args.verbose
    setup_logging(verbosity)

    # Route to appropriate command handler
    if args.command == 'run':
        from qcmanybody.cli.commands.run import handle_run
        return handle_run(args)
    elif args.command == 'plan':
        from qcmanybody.cli.commands.plan import handle_plan
        return handle_plan(args)
    elif args.command == 'validate':
        from qcmanybody.cli.commands.validate import handle_validate
        return handle_validate(args)
    elif args.command == 'convert':
        from qcmanybody.cli.commands.convert import handle_convert
        return handle_convert(args)

if __name__ == '__main__':
    sys.exit(main())
```

**Note**: Using argparse from the standard library keeps dependencies minimal while providing all necessary CLI functionality.

### 2. Input Parser (`input_parser.py`)

**Input File Schema (YAML):**

```yaml
# Basic example
molecule:
  source: xyz  # or inline, qcschema, pdb
  file: water_trimer.xyz
  # OR for inline:
  # source: inline
  # inline:
  #   symbols: [O, H, H, O, H, H, O, H, H]
  #   geometry: [[...], [...], ...]
  #   fragments: [[0,1,2], [3,4,5], [6,7,8]]
  #   units: angstrom  # or bohr

calculation:
  type: single  # or multi
  single:
    driver: gradient  # energy, gradient, hessian
    method: mp2
    basis: cc-pvdz
    program: psi4

  # OR for multi-level:
  # type: multi
  # multi:
  #   driver: energy
  #   levels:
  #     1: {method: ccsd, basis: cc-pvtz, program: psi4}
  #     2: {method: mp2, basis: cc-pvdz, program: psi4}

bsse:
  type: [cp]  # nocp, cp, vmfc - can be list: [cp, vmfc]

manybody:
  max_nbody: 3  # Maximum n-body level
  return_total_data: true
  supersystem_ie_only: false

output:
  file: results.json
  format: json  # json, yaml, text
```

**Parser Functions:**
```python
def parse_input_file(filepath: str) -> dict:
    """Parse YAML/JSON input file"""

def validate_input_schema(data: dict) -> None:
    """Validate against schema, raise clear errors"""

def convert_to_manybody_input(data: dict) -> ManyBodyInput:
    """Convert to ManyBodyInput pydantic model"""

def load_molecule(mol_spec: dict) -> Molecule:
    """Load molecule from various sources"""
```

### 3. Output Writer (`output_writer.py`)

**Responsibilities:**
- Format output in requested format
- Generate summary reports
- Write results to file or stdout

**Output Formats:**

```python
class OutputFormatter(ABC):
    @abstractmethod
    def format(self, result: ManyBodyResult) -> str:
        pass

class JSONFormatter(OutputFormatter):
    """JSON output (machine-readable)"""

class YAMLFormatter(OutputFormatter):
    """YAML output (human-readable)"""

class TextFormatter(OutputFormatter):
    """Human-readable text summary"""

class HDF5Formatter(OutputFormatter):
    """HDF5 format for large data (future)"""
```

## Command-Line Argument Structure

### Global Options
```bash
--verbose, -v           # Increase verbosity (repeatable)
--quiet, -q            # Suppress non-essential output
--version              # Show version
--help, -h             # Show help
--config FILE          # Load default options from file
```

### `run` Command
```bash
qcmanybody run INPUT [OPTIONS]

Required:
  INPUT                 # Path to input file (.yaml or .json)

Options:
  -o, --output FILE     # Output file path (default: stdout)
  --format FORMAT       # Output format: json, yaml, text (default: json)
  --mode MODE           # Execution mode: run, debug, dry-run (default: run)

Logging:
  --log FILE            # Log file path
  --log-level LEVEL     # Log level: debug, info, warning, error
```

### `plan` Command
```bash
qcmanybody plan INPUT [OPTIONS]

Shows execution plan without running calculations.

Options:
  --show-tasks          # List all tasks to be computed
```

### `validate` Command
```bash
qcmanybody validate INPUT

Validates input file without execution.

Options:
  --strict              # Enable strict validation
  --show-schema         # Show expected schema
```

### `convert` Command
```bash
qcmanybody convert INPUT OUTPUT [OPTIONS]

Convert between input formats or from Python script.

Options:
  --from FORMAT         # Source format: yaml, json, python
  --to FORMAT           # Target format: yaml, json
```

## Error Handling Strategy

### Input Validation Errors
- Parse errors: Show line number and context
- Schema errors: Show expected vs. actual with helpful message
- File not found: Show search paths and suggestions

### Execution Errors
- Task failures: Log details, continue with other tasks if possible
- Timeout: Save checkpoint, report progress
- Resource errors: Clear message about limits

### Recovery Mechanisms
- Automatic retry for transient failures
- Checkpoint/resume for long calculations
- Partial result saving

## Integration with pyproject.toml

```toml
[project.scripts]
qcmanybody = "qcmanybody.cli.main:main"

[project.optional-dependencies]
cli = [
    "pyyaml>=5.0",  # Optional: for YAML input file support (JSON works without it)
    "rich>=10.0",   # Optional: for enhanced terminal output (improves UX)
]
```

**Notes:**
- No required dependencies for CLI - uses stdlib `argparse`, `json`, `logging`, `pathlib`
- PyYAML is optional - CLI works with JSON-only if not installed
- Rich is optional - provides prettier output but not required
- Uses existing QCManyBody dependencies: `pydantic`, `qcelemental`, `qcengine`

## Testing Strategy

### Unit Tests
- Test each command in isolation
- Test input parsing with various formats
- Test output formatting
- Test error handling

### Integration Tests
- Full end-to-end CLI execution
- Test with example input files

### Performance Tests
- Compare CLI vs. direct API performance
- Memory usage profiling

## Future Enhancements

1. **Interactive Mode**: TUI for interactive calculation setup
2. **Job Submission**: Integration with SLURM, PBS, SGE
3. **Result Database**: Store and query historical results
4. **Configuration Profiles**: Named configurations for common setups
5. **Shell Completion**: Bash/Zsh completion scripts
6. **Progress Bars**: Rich progress indicators for long calculations
