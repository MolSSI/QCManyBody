# QCManyBody CLI Addition - Project Overview

## Executive Summary

This project adds a comprehensive command-line interface (CLI) to QCManyBody, enabling users to run many-body expansion calculations directly from the terminal using input files, rather than writing Python scripts.

## Current State

QCManyBody currently operates as a Python library with two main interfaces:
- **ManyBodyCore**: Low-level interface for processing pre-computed results
- **ManyBodyComputer**: High-level interface that orchestrates QC calculations via QCEngine

Users must write Python scripts to:
1. Create molecular systems
2. Configure calculation parameters
3. Run calculations
4. Process results

## Target State

The CLI will enable users to:
```bash
# Basic usage
qcmanybody run input.json -o output.json

# With verbose output
qcmanybody run input.json --verbose

# Dry-run mode
qcmanybody plan input.json --show-tasks

# Validate input before running
qcmanybody validate input.json
```

## Key Features

### 1. Input File Support
- JSON input format (primary, stdlib only)
- YAML input format (optional, requires PyYAML)
- Support for ManyBodyComputer workflow via ManyBodyInput
- Validation with clear error messages

### 2. Command Structure
- `qcmanybody run` - Execute calculations
- `qcmanybody plan` - Show calculation plan without running
- `qcmanybody validate` - Validate input file
- `qcmanybody convert` - Convert between JSON and YAML formats

### 3. Output Options
- JSON output (default)
- YAML output
- Text summary output
- Configurable output destination (file or stdout)

### 4. Common Flags
- Input/output file paths
- Logging verbosity levels (`-v`, `-vv`, `-q`)
- Output format selection
- Help and version information

## Technical Approach

### Implementation Strategy
1. Use Python's `argparse` or modern alternative (`click` or `typer`)
2. Integrate with existing ManyBodyComputer and ManyBodyCore
3. Minimal changes to core library code
4. CLI code lives in `qcmanybody/cli/` directory
5. Entry point registered in `pyproject.toml`

### Input File Format
JSON/YAML configuration that maps to ManyBodyInput schema:
- Molecule specification (inline, XYZ, PDB, or QCSchema JSON file)
- Method/basis specifications (single-level or multi-level)
- BSSE treatment options
- Driver (energy, gradient, hessian)
- Many-body expansion options
- Output preferences

### Output Formats
- JSON (default, machine-readable)
- YAML (human-readable)
- Text summary format

## Success Criteria

1. Users can run calculations from command line without writing Python
2. Input validation provides clear, actionable error messages
3. CLI integrates seamlessly with existing ManyBodyComputer
4. Zero new required dependencies (uses only Python stdlib)
5. Comprehensive documentation and examples
6. Full test coverage of CLI functionality

## Project Scope

### In Scope
- Core CLI framework and command structure
- Input file parsing and validation (JSON/YAML)
- Integration with ManyBodyComputer
- Output formatting and export
- Molecule loading from multiple formats
- Documentation and examples
- Test suite

### Out of Scope (Separate Development)
- **Parallel execution** - being developed in separate branch
- GUI or web interface
- Interactive TUI (text user interface)
- Cloud/HPC job submission integration
- Advanced workflow management
- Result database/storage system
- Checkpoint/resume functionality

## Dependencies

### Required
- **No new required dependencies**
- CLI uses only Python standard library (`argparse`, `json`, `logging`, `pathlib`)
- Uses existing QCManyBody dependencies: `pydantic`, `qcelemental`, `numpy`
- Uses existing optional dependency: `qcengine` (for high-level interface)

### Optional
- `pyyaml` - For YAML input file support (JSON works without it)
- `rich` - For enhanced terminal output (improves UX, not yet implemented)

**Philosophy**: The CLI works out-of-the-box with zero additional dependencies beyond what QCManyBody already requires. Optional packages enhance functionality but are not required for core features.

## Timeline Estimate

- **Phase 1**: Planning and Design (Current)
- **Phase 2**: Core CLI Framework (1-2 days)
- **Phase 3**: Input Parsing and Validation (1-2 days)
- **Phase 4**: Integration with ManyBodyComputer (1 day)
- **Phase 5**: Parallel Execution Support (1-2 days)
- **Phase 6**: Testing and Documentation (1-2 days)
- **Total**: ~1 week of development time

## Risks and Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Breaking existing API | High | CLI is additive, no changes to core |
| Input format complexity | Medium | Start simple, iterate based on feedback |
| Parallel execution bugs | Medium | Comprehensive testing, start with simple modes |
| Documentation lag | Low | Document as we build |

## References

- QCManyBody paper: J. Chem. Phys. 161, 152501 (2024) https://doi.org/10.1063/5.0231843
- QCElemental: https://github.com/MolSSI/QCElemental
- QCEngine: https://github.com/MolSSI/QCEngine
