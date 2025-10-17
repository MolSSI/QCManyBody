# QCManyBody CLI - Design Decisions and Rationale

## Overview

This document records key design decisions made during the CLI development, along with rationale and alternatives considered.

---

## Decision 1: CLI Framework

**Decision**: Use `argparse` as the CLI framework

**Status**: ✅ APPROVED

**Date**: 2025-10-15

### Options Considered

1. **argparse** (stdlib) ✅ CHOSEN
   - ✅ No external dependencies
   - ✅ Standard library, stable, well-documented
   - ✅ Fully capable of handling all CLI requirements
   - ✅ Widely understood by Python developers
   - ✅ Battle-tested across thousands of projects
   - ⚠️ More verbose than decorator-based alternatives

2. **click** (third-party)
   - ✅ Clean, decorator-based API
   - ✅ Excellent documentation
   - ✅ Rich features (colors, progress bars, etc.)
   - ✅ Widely used in scientific Python tools
   - ❌ External dependency (not necessary for requirements)
   - ✅ Mature and stable

3. **typer** (third-party)
   - ✅ Modern, type-hint based API
   - ✅ Built on click
   - ✅ Automatic help generation
   - ❌ Newer, less established
   - ❌ Additional dependency

### Rationale

**Decision: argparse**

- **Zero dependencies**: Part of Python standard library, no installation overhead
- **Fully sufficient**: Capable of handling all CLI requirements for this project
- **Standardized**: Well-established Python standard, excellent documentation
- **Maintainability**: Will be supported as long as Python exists
- **Project philosophy**: Avoid unnecessary dependencies when built-in tools are adequate
- **Simplicity**: While more verbose, the code is straightforward and explicit

While click offers a cleaner API, argparse is more than adequate for QCManyBody's needs and avoids adding an external dependency for functionality that's available in the standard library.

**Implementation Example:**

```python
import argparse

def main():
    parser = argparse.ArgumentParser(
        prog='qcmanybody',
        description='QCManyBody: Many-body expansion calculations'
    )

    # Add version
    parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

    # Add subparsers for commands
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Run command
    run_parser = subparsers.add_parser('run', help='Run a calculation')
    run_parser.add_argument('input', type=str, help='Input file path')
    run_parser.add_argument('-o', '--output', help='Output file path')
    run_parser.add_argument('--format',
                          choices=['json', 'yaml', 'text'],
                          default='json',
                          help='Output format')

    # Parse arguments
    args = parser.parse_args()

    # Route to appropriate handler
    if args.command == 'run':
        handle_run(args)
```

This approach is clear, explicit, and requires no external dependencies.

---

## Decision 2: Input File Format

**Decision**: Support JSON (primary, required) and YAML (optional, when PyYAML available)

**Status**: ✅ APPROVED

**Date**: 2025-10-15

### Options Considered

1. **JSON only**
   - ✅ No external dependencies (stdlib json)
   - ✅ Language-agnostic
   - ✅ Easy to generate programmatically
   - ✅ Sufficient for all use cases
   - ❌ No comments
   - ⚠️ Less human-readable than YAML

2. **YAML only**
   - ✅ Human-readable, great for hand-editing
   - ✅ Supports comments
   - ✅ Less verbose than JSON
   - ❌ Requires PyYAML dependency
   - ❌ More complex parsing

3. **Both JSON and YAML (with JSON as primary)** ✅ CHOSEN
   - ✅ JSON always works (no dependencies)
   - ✅ YAML for human editing, JSON for automation
   - ✅ Easy conversion between formats
   - ❌ Two parsers to maintain
   - ❌ PyYAML dependency

4. **TOML**
   - ✅ Modern, clean syntax
   - ✅ Growing popularity in Python
   - ❌ Less familiar to most users
   - ❌ Additional dependency

5. **Custom INI-style**
   - ✅ Simple syntax
   - ❌ Limited expressiveness
   - ❌ Non-standard, poor tooling support

### Rationale

**Decision: JSON primary, YAML optional**

- **JSON (primary, no dependencies)**:
  - Part of Python standard library
  - Sufficient for all use cases
  - Universal format with excellent tooling
  - Easy to generate programmatically
  - Language-agnostic, good for automation
  - Works out of the box with no additional installation

- **YAML (optional enhancement)**:
  - More human-readable for manual editing
  - Supports comments for documentation
  - Less verbose than JSON
  - Common in scientific computing workflows
  - **Only enabled if PyYAML is installed** (graceful fallback to JSON-only)

**Implementation Strategy:**
```python
# Try to import yaml, fall back gracefully
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False
    # CLI still works, just JSON-only

def load_input_file(filepath):
    """Load input file, supporting JSON and YAML if available."""
    if filepath.endswith('.yaml') or filepath.endswith('.yml'):
        if not YAML_AVAILABLE:
            raise RuntimeError(
                "YAML file provided but PyYAML not installed. "
                "Install with: pip install pyyaml\n"
                "Or use JSON format instead."
            )
        with open(filepath) as f:
            return yaml.safe_load(f)
    else:  # JSON
        with open(filepath) as f:
            return json.load(f)
```

**User Experience:**
- Users can start immediately with JSON (no extra dependencies)
- Power users can optionally install PyYAML for YAML support
- Clear error message if YAML file provided without PyYAML
- Documentation shows both formats with JSON as default

**Dependency approach:**
- PyYAML listed in `[project.optional-dependencies]` CLI section
- Not required for basic CLI functionality
- Users choose based on their preference

---

## Decision 3: Schema Validation

**Decision**: Use Pydantic models for input validation

**Status**: Approved

**Date**: 2025-10-15

### Options Considered

1. **Pydantic**
   - ✅ Already a dependency of QCManyBody
   - ✅ Rich validation features
   - ✅ Clear error messages
   - ✅ Type hints
   - ✅ JSON Schema generation

2. **jsonschema**
   - ✅ Standard JSON Schema validation
   - ❌ Additional dependency
   - ❌ Less Pythonic API

3. **Manual validation**
   - ✅ No dependencies
   - ❌ Error-prone
   - ❌ Time-consuming to maintain

### Rationale

**Decision: Pydantic**

- Already a core dependency (no new dependencies)
- Excellent validation capabilities
- Clear, helpful error messages
- Type safety throughout the codebase
- Easy to extend and maintain
- Can generate JSON Schema for documentation

---

## Decision 4: Output Format

**Decision**: Support JSON, YAML, and text summary formats

**Status**: Recommended

**Date**: 2025-10-15

### Options Considered

1. **JSON only**
   - ✅ Machine-readable
   - ✅ Standard format
   - ❌ Not human-friendly

2. **YAML only**
   - ✅ Human-readable
   - ❌ Less standard for output

3. **HDF5**
   - ✅ Efficient for large datasets
   - ✅ Standard in scientific computing
   - ❌ Binary format
   - ❌ Requires h5py dependency
   - ❌ Overkill for most use cases

4. **Multiple formats**
   - ✅ Flexibility for different use cases
   - ✅ JSON for automation, YAML for readability, text for summaries
   - ❌ More code to maintain

### Rationale

**Recommendation: JSON (default), YAML, and text summary**

- **JSON**: Default format for:
  - Automation and scripting
  - Integration with other tools
  - Machine parsing

- **YAML**: Optional format for:
  - Human review of results
  - Documentation and reports

- **Text summary**: Optional format for:
  - Quick inspection of results
  - Publication-ready output
  - Terminal display

- **HDF5**: Defer to future enhancement
  - Add when needed for very large datasets
  - Not needed for typical use cases

---

## Decision 5: Parallel Execution

**Decision**: Deferred to separate branch

**Status**: ❌ NOT IMPLEMENTED IN THIS VERSION

**Date**: 2025-10-16

### Rationale

Parallel execution (multiprocessing, threading, MPI) is being developed in a separate branch and is not part of this CLI implementation. This decision was made to:

- **Minimize scope**: Focus on core CLI functionality first
- **Reduce complexity**: Avoid parallel execution bugs and edge cases
- **Avoid core changes**: Parallel execution may require modifications to existing core code
- **Enable focused development**: Parallel features can be developed and tested independently

The CLI implementation will use QCManyBody's existing serial execution mode via the standard ManyBodyComputer interface.

### Future Work

Parallel execution support will be added in a future version, developed in a separate branch with:
- Multiprocessing executor for workstations
- MPI executor for HPC clusters
- Proper task distribution and result collection
- Integration with the CLI

---

## Decision 6: Error Handling Philosophy

**Decision**: Fail-fast for input errors, continue-on-error for computation

**Status**: Approved

**Date**: 2025-10-15

### Principles

1. **Input validation**: Fail immediately with clear error messages
   - Invalid file format → Error before execution
   - Invalid schema → Error with specific field and suggestion
   - Missing files → Error with helpful message

2. **Computation errors**: Continue with other tasks when possible
   - Task failure → Log error, continue with remaining tasks
   - Partial results → Save successfully computed results
   - Summary report → Show which tasks succeeded/failed

3. **Resource errors**: Fail gracefully with cleanup
   - Out of memory → Clean up, report error, exit
   - Disk space → Clean up temporary files, exit
   - Program errors → Log details, provide diagnostic information

### Rationale

- **User time is valuable**: Don't waste hours of computation due to a single task failure
- **Clear feedback**: Users should know immediately if their input is wrong
- **Robustness**: Partial results are still useful even if some tasks fail
- **Debugging**: Failed tasks are logged with enough detail to diagnose

---

## Decision 7: Progress Reporting

**Decision**: Use structured logging with optional rich terminal output

**Status**: Recommended

**Date**: 2025-10-15

### Options Considered

1. **Print statements**
   - ❌ Unstructured, hard to control
   - ❌ Doesn't integrate with logging systems

2. **Python logging**
   - ✅ Standard, structured
   - ✅ Configurable levels
   - ✅ Can write to file
   - ❌ Plain text only

3. **Rich library**
   - ✅ Beautiful terminal output
   - ✅ Progress bars, tables, colors
   - ✅ Better UX
   - ❌ Additional dependency
   - ❌ Doesn't work well in log files

4. **Logging + optional Rich**
   - ✅ Best of both worlds
   - ✅ Rich for terminal, logging for files
   - ✅ Graceful fallback

### Rationale

**Recommendation: Python logging with optional Rich enhancement**

- Use standard logging for all messages
- Optionally enhance terminal output with Rich when available
- Automatically detect if output is a TTY
- In non-TTY (log files, CI), use plain text

Example:
```python
# Core: always available
logger.info("Task 5/24 completed")

# Enhanced: only if Rich available and stdout is TTY
if rich_available and sys.stdout.isatty():
    progress.update(task_id, completed=5, total=24)
```

---

## Decision 8: Configuration File Support

**Decision**: Support optional configuration files in YAML/JSON format

**Status**: Recommended (defer to future version)

**Date**: 2025-10-15

### Use Cases

Users may want default settings:
- Default program (psi4 vs nwchem)
- Default output format
- Default verbosity level

### Configuration File Locations

Search order:
1. `--config` flag (explicit)
2. `./qcmanybody.yaml` (project-specific)
3. `~/.config/qcmanybody/config.yaml` (user-specific)
4. System defaults

### Rationale

- Useful for power users and HPC environments
- Command-line flags override config file
- Input file overrides config file
- Precedence: CLI flags > Input file > Config file > Defaults

---

## Decision 9: Checkpointing and Resume

**Decision**: Deferred to separate branch

**Status**: ❌ NOT IMPLEMENTED IN THIS VERSION

**Date**: 2025-10-16

### Rationale

Checkpointing and resume functionality is being developed in a separate branch along with parallel execution features. This decision was made to:

- **Minimize scope**: Focus on core CLI functionality
- **Reduce complexity**: Checkpointing adds significant state management complexity
- **Avoid dependencies**: Different checkpoint formats have trade-offs
- **Enable focused development**: Can be developed with parallel execution features

### Future Work

Checkpoint/resume support will be added in a future version with:
- Periodic checkpoint saves during calculation
- Resume from checkpoint file
- Atomic writes for data integrity
- Version compatibility checking

---

## Decision 10: Minimal Core Modifications

**Decision**: CLI is purely additive, no changes to core library

**Status**: Approved

**Date**: 2025-10-15

### Principles

1. All CLI code lives in `qcmanybody/cli/` directory
2. No modifications to existing APIs
3. CLI uses public interfaces only
4. Core library can be used without CLI

### Rationale

- Reduces risk of breaking existing users
- Clear separation of concerns
- CLI can evolve independently
- Easier to review and test

### Exception

If minor utility functions would benefit both CLI and core, add them to a shared `utils.py` with careful consideration.

---

## Decision 11: Testing Strategy

**Decision**: Comprehensive testing at multiple levels

**Status**: Approved

**Date**: 2025-10-15

### Test Levels

1. **Unit tests**: Test individual functions and classes
   - Input parser functions
   - Output formatters
   - Executor logic

2. **Integration tests**: Test end-to-end flows
   - Full calculation runs
   - Different input formats
   - Error handling paths

3. **CLI tests**: Test command-line interface
   - Argument parsing
   - Help text
   - Exit codes

4. **Regression tests**: Compare to Python API
   - Results must match exactly
   - Performance must be comparable

### Test Data

- Use existing test fixtures from QCManyBody
- Create small, fast test cases for CLI-specific features
- Use mocking for QCEngine calls when appropriate

### Rationale

- High test coverage ensures quality
- Multiple test levels catch different issues
- Regression tests ensure CLI matches Python API exactly

---

## Decision 12: Documentation Structure

**Decision**: Integrate CLI docs into existing MkDocs site

**Status**: Recommended

**Date**: 2025-10-15

### Documentation Components

1. **User Guide**: Getting started, basic usage
2. **Reference**: Command-line options, input file format
3. **Examples**: Common use cases
4. **Migration Guide**: From Python scripts to CLI
5. **API Docs**: For CLI module (auto-generated)

### Integration

- Add "CLI Usage" section to main docs
- Keep existing Python API docs unchanged
- Cross-reference between CLI and API docs

### Rationale

- Centralized documentation
- Easy to maintain
- Searchable
- Version-controlled

---

## Decision 13: Deprecation and Versioning

**Decision**: CLI follows QCManyBody version, separate feature versioning via schema

**Status**: Approved

**Date**: 2025-10-15

### Version Strategy

- CLI version follows QCManyBody package version
- Input file schema has separate version (starting at 1.0)
- CLI must handle multiple schema versions (backward compatibility)

### Deprecation Policy

- Features deprecated in minor versions
- Removed in major versions
- Clear warning messages
- Migration guides provided

### Rationale

- Users don't need to track separate CLI version
- Schema versioning allows input file evolution
- Backward compatibility prevents breakage

---

## Future Considerations

### Items for Future Discussion

1. **Interactive TUI Mode**
   - Would users benefit from an interactive mode?
   - Consider `textual` or `urwid` libraries

2. **Job Submission Integration**
   - Direct SLURM/PBS/SGE submission from CLI?
   - Or keep as separate user responsibility?

3. **Result Database**
   - Should CLI support storing results in a database?
   - SQLite, PostgreSQL, or dedicated scientific DB?

4. **Cloud Execution**
   - Support for cloud platforms (AWS, GCP, Azure)?
   - Likely out of scope for v1

5. **Web Interface**
   - Complement CLI with web UI?
   - Separate project or integrated?

6. **Plugin System**
   - Allow third-party extensions?
   - Custom executors, formatters, validators?

---

## Review and Approval

| Decision | Reviewer | Status | Date | Notes |
|----------|----------|--------|------|-------|
| CLI Framework (argparse) | User | ✅ Approved | 2025-10-15 | Use stdlib, avoid unnecessary deps |
| Input Format (JSON primary, YAML optional) | User | ✅ Approved | 2025-10-15 | JSON required, YAML optional |
| Schema Validation (Pydantic) | - | ✅ Approved | 2025-10-15 | Already a dependency |
| Output Formats | - | ✅ Approved | 2025-10-15 | JSON/YAML/text |
| Parallel Execution | User | ❌ Deferred | 2025-10-16 | Separate branch, not in this version |
| Error Handling | - | ✅ Approved | 2025-10-15 | Fail-fast input, continue-on-error compute |
| Progress Reporting | - | ✅ Approved | 2025-10-15 | Logging with optional rich |
| Configuration Files | - | Deferred | - | Future version |
| Checkpointing | User | ❌ Deferred | 2025-10-16 | Separate branch with parallel execution |
| Minimal Core Modifications | User | ✅ Approved | 2025-10-16 | CLI is purely additive |
| Testing Strategy | - | ✅ Approved | 2025-10-15 | Multi-level comprehensive testing |

---

**Last Updated**: 2025-10-16
**Next Review**: During implementation
