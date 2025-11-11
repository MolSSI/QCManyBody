# QCManyBody CLI Addition - Project Documentation

## Welcome

This directory contains all planning, design, and tracking documentation for the QCManyBody command-line interface (CLI) addition project.

## Project Goal

Add a comprehensive command-line interface to QCManyBody, enabling users to run many-body expansion calculations directly from the terminal using input files, rather than writing Python scripts.

## Documentation Structure

### Planning Documents

1. **[OVERVIEW.md](OVERVIEW.md)** - Project overview and high-level goals
   - Executive summary
   - Current vs. target state
   - Key features
   - Success criteria
   - Timeline and scope

2. **[ARCHITECTURE.md](ARCHITECTURE.md)** - Technical design and architecture
   - System architecture diagrams
   - Module structure
   - Component responsibilities
   - Implementation details
   - Command-line argument structure

3. **[TASKS.md](TASKS.md)** - Detailed task breakdown
   - Phase-by-phase task lists
   - Task dependencies
   - Milestone checklists
   - Future enhancements

4. **[PROGRESS.md](PROGRESS.md)** - Project tracking and status
   - Phase progress tracking
   - Recent activity log
   - Sprint planning
   - Milestone tracking
   - Risk and issue management
   - Test coverage metrics

### Specification Documents

5. **[INPUT_FILE_SPEC.md](INPUT_FILE_SPEC.md)** - Input file format specification
   - Complete schema documentation
   - Field descriptions and options
   - Validation rules
   - Examples for each section
   - Best practices

6. **[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md)** - Key design decisions and rationale
   - Decision log with rationale
   - Alternatives considered
   - Trade-off analysis
   - Approval status

### Usage Documentation

7. **[EXAMPLES.md](EXAMPLES.md)** - Usage examples and tutorials
   - Quick start guide
   - Basic usage examples
   - Advanced features
   - HPC integration
   - Batch processing
   - Error handling

## Quick Reference

### Current Status

- **Phase**: Phase 1 - Foundation and Planning
- **Progress**: ~80% of planning complete
- **Status**: Planning documents created, awaiting review

### Key Decisions Pending

1. **CLI Framework**: argparse vs. click (Recommendation: click)
2. **Input Format**: YAML, JSON, or both (Recommendation: both)
3. **Optional Dependencies**: Which to include (Recommendation: click + pyyaml)

### Next Steps

1. Review planning documents
2. Finalize key design decisions
3. Begin Phase 2: Core CLI Framework implementation
4. Update pyproject.toml with CLI dependencies

## Project Structure

Planned implementation structure:

```
qcmanybody/
├── cli/                          # New CLI package
│   ├── __init__.py
│   ├── main.py                   # Entry point
│   ├── commands/                 # Command implementations
│   │   ├── __init__.py
│   │   ├── run.py
│   │   ├── plan.py
│   │   ├── validate.py
│   │   └── convert.py
│   ├── input_parser.py           # Input file parsing
│   ├── executor.py               # Execution coordination
│   ├── output_writer.py          # Output formatting
│   ├── schemas/                  # Input schemas and examples
│   │   ├── __init__.py
│   │   ├── input_schema.py
│   │   └── examples/
│   └── utils.py
├── core.py                        # Existing (unchanged)
├── computer.py                    # Existing (unchanged)
└── ...
```

## How to Use This Documentation

### For Project Planning
1. Start with **OVERVIEW.md** for the big picture
2. Review **TASKS.md** for detailed work breakdown
3. Track progress in **PROGRESS.md**

### For Implementation
1. Reference **ARCHITECTURE.md** for technical design
2. Use **INPUT_FILE_SPEC.md** for input parsing implementation
3. Consult **DESIGN_DECISIONS.md** for rationale behind choices

### For Testing
1. Use **EXAMPLES.md** for test cases
2. Refer to **INPUT_FILE_SPEC.md** for validation test scenarios

### For Documentation
1. **EXAMPLES.md** provides user-facing usage examples
2. **INPUT_FILE_SPEC.md** becomes reference documentation
3. **OVERVIEW.md** can be adapted for user guide introduction

## Development Phases

### ✅ Phase 1: Foundation and Planning (Current)
- Project setup
- Planning documentation
- Design decisions
- Example creation

### ⏳ Phase 2: Core CLI Framework (Next)
- Package structure
- Entry point implementation
- Command stubs
- Basic functionality

### ⏳ Phase 3: Input Parsing
- Schema definition
- Parser implementation
- Validation
- Error messages

### ⏳ Phase 4: ManyBodyComputer Integration
- Executor implementation
- Run command
- Plan command
- Validate command

### ⏳ Phase 5: Output Formatting
- Format implementations
- Output writers
- Summary generation

### ⏳ Phase 6: Parallel Execution
- Executor variants
- Performance testing
- Scaling validation

### ⏳ Phase 7: Advanced Features
- Checkpointing
- Configuration files
- Additional utilities

### ⏳ Phase 8: Testing
- Unit tests
- Integration tests
- CLI tests
- Performance tests

### ⏳ Phase 9: Documentation
- User guide
- API docs
- Tutorials
- Migration guide

### ⏳ Phase 10: QA and Release
- Code quality
- Performance benchmarking
- User testing
- Release preparation

## Contributing to This Documentation

### Adding New Documents
- Place in `cli_addition/` directory
- Update this README with link and description
- Follow existing formatting style

### Updating Existing Documents
- Keep version history at bottom of document
- Note date and nature of changes
- Update "Last Updated" timestamp

### Document Templates

#### Planning Document Template
```markdown
# Title

## Overview
[Brief description]

## Details
[Main content]

## References
[Links to related docs]

---
**Last Updated**: YYYY-MM-DD
```

## Resources

### External References
- [QCManyBody GitHub](https://github.com/MolSSI/QCManyBody)
- [QCManyBody Docs](https://molssi.github.io/QCManyBody/)
- [QCManyBody Paper](https://doi.org/10.1063/5.0231843) - J. Chem. Phys. 161, 152501 (2024)
- [QCElemental](https://github.com/MolSSI/QCElemental)
- [QCEngine](https://github.com/MolSSI/QCEngine)

### Related Tools
- [click](https://click.palletsprojects.com/) - CLI framework
- [Pydantic](https://docs.pydantic.dev/) - Validation
- [Rich](https://rich.readthedocs.io/) - Terminal output

## Questions and Feedback

For questions or feedback about the CLI design:
1. Review relevant planning documents
2. Check **DESIGN_DECISIONS.md** for rationale
3. Open a discussion with project maintainers

## License

This planning documentation follows the same license as QCManyBody (BSD-3-Clause).

---

**Project Start**: 2025-10-15
**Current Phase**: Phase 1 - Planning
**Last Updated**: 2025-10-15
