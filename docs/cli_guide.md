# QCManyBody CLI User Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Commands](#commands)
5. [Input File Format](#input-file-format)
6. [Examples](#examples)
7. [Troubleshooting](#troubleshooting)

---

## Introduction

The QCManyBody command-line interface (CLI) provides a user-friendly way to run many-body expansion calculations without writing Python code. The CLI accepts input files in JSON or YAML format and executes calculations through the QCManyBody library.

### Key Features

- **Simple input format**: Human-readable JSON or YAML files
- **Four core commands**: `run`, `plan`, `validate`, `convert`
- **Multiple output formats**: JSON, YAML, or text summary
- **No code required**: Perfect for scripting and batch processing
- **Comprehensive validation**: Catch errors before computation starts

---

## Installation

### Install QCManyBody

```bash
# Install from PyPI (when released)
pip install qcmanybody[cli]

# Or install from source with CLI support
git clone https://github.com/MolSSI/QCManyBody.git
cd QCManyBody
pip install -e .
```

### Optional Dependencies

```bash
# For YAML support (highly recommended)
pip install pyyaml

# For QC program execution
pip install qcengine

# Specific QC programs (at least one required for calculations)
conda install -c psi4 psi4
# or
conda install -c conda-forge nwchem
```

### Verify Installation

```bash
qcmanybody --version
qcmanybody --help
```

---

## Quick Start

### 1. Create an Input File

Create a file called `water_dimer.json`:

```json
{
  "schema_name": "qcmanybody_cli_input",
  "schema_version": 1,
  "molecule": {
    "source": "inline",
    "inline": {
      "symbols": ["O", "H", "H", "O", "H", "H"],
      "geometry": [
        [0.000, 0.000, 0.000],
        [0.757, 0.586, 0.000],
        [-0.757, 0.586, 0.000],
        [3.000, 0.000, 0.000],
        [3.757, 0.586, 0.000],
        [2.243, 0.586, 0.000]
      ],
      "fragments": [[0, 1, 2], [3, 4, 5]],
      "units": "angstrom"
    }
  },
  "calculation": {
    "type": "single",
    "single": {
      "driver": "energy",
      "method": "scf",
      "basis": "sto-3g",
      "program": "psi4"
    }
  },
  "bsse": {
    "type": ["cp"]
  },
  "manybody": {
    "max_nbody": 2
  }
}
```

### 2. Validate Your Input

```bash
qcmanybody validate water_dimer.json
```

### 3. Preview the Calculation Plan

```bash
qcmanybody plan water_dimer.json
```

### 4. Run the Calculation

```bash
qcmanybody run water_dimer.json -o results.json
```

---

## Commands

### `qcmanybody run`

Execute a many-body calculation.

**Usage:**
```bash
qcmanybody run <input-file> [options]
```

**Options:**
- `-o, --output FILE`: Write output to file (default: stdout)
- `--format {json,yaml,text}`: Output format (default: json)
- `--log FILE`: Write log messages to file

**Examples:**
```bash
# Run with JSON output
qcmanybody run input.json -o results.json

# Run with text summary
qcmanybody run input.yaml --format text

# Run with logging
qcmanybody run input.json --log calculation.log
```

---

### `qcmanybody plan`

Show what calculations will be performed without executing them.

**Usage:**
```bash
qcmanybody plan <input-file> [options]
```

**Options:**
- `--show-tasks`: Display detailed task breakdown
- `--log FILE`: Write log messages to file

**Examples:**
```bash
# Show basic execution plan
qcmanybody plan input.json

# Show detailed task list
qcmanybody plan input.json --show-tasks
```

**Output includes:**
- Molecular system information (atoms, fragments)
- Calculation settings (method, basis, program)
- BSSE correction types
- Total number of tasks
- Task breakdown by n-body level (with --show-tasks)

---

### `qcmanybody validate`

Validate an input file without running calculations.

**Usage:**
```bash
qcmanybody validate <input-file> [options]
qcmanybody validate --show-schema
```

**Options:**
- `--show-schema`: Display the JSON schema for input files
- `--strict`: Enable strict validation mode
- `--log FILE`: Write log messages to file

**Examples:**
```bash
# Validate an input file
qcmanybody validate input.json

# Show expected schema
qcmanybody validate --show-schema
```

**Validation checks:**
- Schema compliance
- Molecule specification correctness
- Fragment definitions (no gaps, no overlaps)
- Calculation settings validity
- BSSE type compatibility
- Conversion to internal format

---

### `qcmanybody convert`

Convert input files between JSON and YAML formats.

**Usage:**
```bash
qcmanybody convert <input-file> <output-file> [options]
```

**Options:**
- `--log FILE`: Write log messages to file

**Examples:**
```bash
# Convert JSON to YAML
qcmanybody convert input.json input.yaml

# Convert YAML to JSON
qcmanybody convert input.yaml input.json
```

**Note:** Output format is determined by the output file extension (`.json`, `.yaml`, or `.yml`).

---

## Input File Format

### Schema Overview

All input files must include:

```json
{
  "schema_name": "qcmanybody_cli_input",
  "schema_version": 1,
  "molecule": { ... },
  "calculation": { ... },
  "bsse": { ... },
  "manybody": { ... }
}
```

### Molecule Specification

#### Inline Molecule

Specify molecular geometry directly:

```json
"molecule": {
  "source": "inline",
  "inline": {
    "symbols": ["He", "He", "He"],
    "geometry": [
      [0.0, 0.0, 0.0],
      [0.0, 0.0, 3.0],
      [0.0, 0.0, 6.0]
    ],
    "fragments": [[0], [1], [2]],
    "molecular_charge": 0.0,
    "molecular_multiplicity": 1,
    "units": "angstrom"
  }
}
```

**Fields:**
- `symbols`: List of element symbols
- `geometry`: List of [x, y, z] coordinates
- `fragments`: List of atom indices for each fragment (0-indexed)
- `molecular_charge`: Total charge (default: 0.0)
- `molecular_multiplicity`: Spin multiplicity 2S+1 (default: 1)
- `units`: "angstrom" or "bohr" (default: "angstrom")
- `fragment_charges`: Charges per fragment (optional)
- `fragment_multiplicities`: Multiplicities per fragment (optional)

#### From File

Load molecule from external file:

```json
"molecule": {
  "source": "file",
  "file": {
    "file": "path/to/molecule.xyz",
    "fragments": [[0, 1, 2], [3, 4, 5]]
  }
}
```

**Supported formats:** XYZ, PDB, QCSchema JSON

---

### Calculation Specification

#### Single-Level Calculation

Use the same method for all n-body levels:

```json
"calculation": {
  "type": "single",
  "single": {
    "driver": "energy",
    "method": "mp2",
    "basis": "cc-pvdz",
    "program": "psi4"
  }
}
```

**Fields:**
- `driver`: "energy", "gradient", or "hessian"
- `method`: QC method (e.g., "scf", "mp2", "ccsd(t)")
- `basis`: Basis set name
- `program`: QC program ("psi4", "nwchem", "cfour")

#### Multi-Level Calculation

Use different methods for different n-body levels:

```json
"calculation": {
  "type": "multi",
  "multi": {
    "driver": "energy",
    "levels": {
      "1": {
        "method": "ccsd(t)",
        "basis": "cc-pvtz",
        "program": "nwchem"
      },
      "2": {
        "method": "mp2",
        "basis": "cc-pvdz",
        "program": "psi4"
      },
      "3": {
        "method": "scf",
        "basis": "cc-pvdz",
        "program": "psi4"
      }
    }
  }
}
```

**Keys:** Level numbers ("1", "2", "3", ...) or "supersystem"

---

### BSSE Correction

Specify basis set superposition error correction:

```json
"bsse": {
  "type": ["cp"]
}
```

**Options:**
- `"nocp"`: No counterpoise correction
- `"cp"`: Counterpoise correction (Boys-Bernardi)
- `"vmfc"`: Valiron-Mayer function counterpoise

**Multiple BSSE types:** The first in the list determines the returned result.

---

### Many-Body Settings

Configure the many-body expansion:

```json
"manybody": {
  "max_nbody": 3,
  "return_total_data": true,
  "supersystem_ie_only": false
}
```

**Fields:**
- `max_nbody`: Maximum n-body level to compute (default: number of fragments)
- `return_total_data`: Return total energy/properties (default: false)
- `supersystem_ie_only`: Only compute supersystem interaction energy (default: false)
- `embedding_charges`: Point charges for embedding (optional)

---

### Program Keywords

Pass program-specific keywords:

```json
"program": {
  "keywords": {
    "scf__e_convergence": 1e-8,
    "mp2__freeze_core": true
  }
}
```

**Format:** Double underscore separates module from keyword (e.g., `scf__e_convergence`)

---

### Output Options

Control output formatting:

```json
"output": {
  "format": "json",
  "pretty": true,
  "include_timings": true
}
```

**Fields:**
- `format`: "json", "yaml", or "text" (default: "json")
- `pretty`: Pretty-print output (default: true)
- `include_timings`: Include timing information (default: true)

---

## Examples

### Example 1: Basic Energy Calculation

Simple He₃ system with CP-corrected energy:

```json
{
  "schema_name": "qcmanybody_cli_input",
  "schema_version": 1,
  "molecule": {
    "source": "inline",
    "inline": {
      "symbols": ["He", "He", "He"],
      "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0], [0.0, 0.0, 6.0]],
      "fragments": [[0], [1], [2]],
      "units": "angstrom"
    }
  },
  "calculation": {
    "type": "single",
    "single": {
      "driver": "energy",
      "method": "mp2",
      "basis": "cc-pvdz",
      "program": "psi4"
    }
  },
  "bsse": {
    "type": ["cp"]
  },
  "manybody": {
    "max_nbody": 3
  }
}
```

**Run:**
```bash
qcmanybody run example1.json -o results1.json
```

---

### Example 2: Gradient Calculation

Water trimer with gradient computation:

```json
{
  "schema_name": "qcmanybody_cli_input",
  "schema_version": 1,
  "molecule": {
    "source": "inline",
    "inline": {
      "symbols": ["O", "H", "H", "O", "H", "H", "O", "H", "H"],
      "geometry": [
        [0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0],
        [3.0, 0.0, 0.0], [3.757, 0.586, 0.0], [2.243, 0.586, 0.0],
        [1.5, 2.598, 0.0], [2.257, 3.184, 0.0], [0.743, 3.184, 0.0]
      ],
      "fragments": [[0, 1, 2], [3, 4, 5], [6, 7, 8]],
      "units": "angstrom"
    }
  },
  "calculation": {
    "type": "single",
    "single": {
      "driver": "gradient",
      "method": "scf",
      "basis": "sto-3g",
      "program": "psi4"
    }
  },
  "manybody": {
    "max_nbody": 2
  }
}
```

**Run:**
```bash
qcmanybody run example2.json --format text
```

---

### Example 3: Multi-Level Calculation

High-accuracy Ne₃ with multi-level expansion:

```json
{
  "schema_name": "qcmanybody_cli_input",
  "schema_version": 1,
  "molecule": {
    "source": "inline",
    "inline": {
      "symbols": ["Ne", "Ne", "Ne"],
      "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 4.0]],
      "fragments": [[0], [1], [2]],
      "units": "angstrom"
    }
  },
  "calculation": {
    "type": "multi",
    "multi": {
      "driver": "energy",
      "levels": {
        "1": {"method": "ccsd(t)", "basis": "cc-pvtz", "program": "nwchem"},
        "2": {"method": "mp2", "basis": "cc-pvdz", "program": "psi4"},
        "3": {"method": "scf", "basis": "cc-pvdz", "program": "psi4"}
      }
    }
  },
  "bsse": {
    "type": ["cp", "vmfc"]
  },
  "manybody": {
    "max_nbody": 3,
    "return_total_data": true
  }
}
```

**Run:**
```bash
qcmanybody plan example3.json --show-tasks
qcmanybody run example3.json -o results3.json
```

---

### Example 4: Loading from XYZ File

Using an external geometry file:

**molecule.xyz:**
```
6
Water dimer
O    0.000    0.000    0.000
H    0.757    0.586    0.000
H   -0.757    0.586    0.000
O    3.000    0.000    0.000
H    3.757    0.586    0.000
H    2.243    0.586    0.000
```

**input.json:**
```json
{
  "schema_name": "qcmanybody_cli_input",
  "schema_version": 1,
  "molecule": {
    "source": "file",
    "file": {
      "file": "molecule.xyz",
      "fragments": [[0, 1, 2], [3, 4, 5]]
    }
  },
  "calculation": {
    "type": "single",
    "single": {
      "driver": "energy",
      "method": "scf",
      "basis": "sto-3g",
      "program": "psi4"
    }
  },
  "bsse": {
    "type": ["cp"]
  },
  "manybody": {
    "max_nbody": 2
  }
}
```

---

## Troubleshooting

### Common Errors

#### 1. "Program X is not available"

**Problem:** QC program not installed or not in PATH.

**Solution:**
```bash
# Check if program is available
which psi4
which nwchem

# Install via conda
conda install -c psi4 psi4
```

#### 2. "Validation failed: field required"

**Problem:** Required field missing from input file.

**Solution:** Run validation to see exactly what's missing:
```bash
qcmanybody validate --show-schema > schema.json
qcmanybody validate input.json
```

#### 3. "Fragment indices out of range"

**Problem:** Fragment definitions reference atoms that don't exist.

**Solution:** Ensure fragment indices are 0-based and within range:
```python
# For 6 atoms, valid indices are 0-5
"fragments": [[0, 1, 2], [3, 4, 5]]  # ✓ Correct
"fragments": [[1, 2, 3], [4, 5, 6]]  # ✗ Wrong (6 is out of range)
```

#### 4. "YAML parsing failed"

**Problem:** PyYAML not installed.

**Solution:**
```bash
pip install pyyaml
```

Or use JSON format instead.

#### 5. "Cannot parse geometry file"

**Problem:** File format not recognized or file not found.

**Solution:**
- Use absolute path or relative path from where you run the command
- Ensure file format is supported (XYZ, PDB, QCSchema JSON)
- Check file permissions

---

### Getting Help

#### Command Help

```bash
# General help
qcmanybody --help

# Command-specific help
qcmanybody run --help
qcmanybody plan --help
qcmanybody validate --help
qcmanybody convert --help
```

#### Show Expected Schema

```bash
qcmanybody validate --show-schema
```

#### Enable Verbose Logging

```bash
qcmanybody -v run input.json
qcmanybody -vv run input.json  # More verbose
qcmanybody -vvv run input.json  # Most verbose
```

---

### Tips and Best Practices

1. **Always validate first:**
   ```bash
   qcmanybody validate input.json && qcmanybody run input.json
   ```

2. **Preview before running:**
   ```bash
   qcmanybody plan input.json --show-tasks
   ```

3. **Use YAML for human readability:**
   - Add comments with `#`
   - More forgiving syntax
   - Convert to JSON for distribution: `qcmanybody convert input.yaml input.json`

4. **Save logs for long calculations:**
   ```bash
   qcmanybody run input.json --log calculation.log
   ```

5. **Start small:**
   - Test with 2-fragment systems first
   - Use small basis sets (sto-3g) for testing
   - Scale up once workflow is validated

6. **Check task count:**
   ```bash
   qcmanybody plan input.json  # Shows total tasks
   ```

   Many-body calculations scale as C(n,k) where n=fragments, k=max_nbody:
   - 3 fragments, max_nbody=2: 6 tasks
   - 4 fragments, max_nbody=3: 14 tasks
   - 5 fragments, max_nbody=4: 30 tasks

---

## Additional Resources

- **GitHub Repository:** https://github.com/MolSSI/QCManyBody
- **Documentation:** https://molssi.github.io/QCManyBody
- **Report Issues:** https://github.com/MolSSI/QCManyBody/issues
- **Example Files:** See `examples/cli/` in the repository

---

## Version Information

This guide is for QCManyBody CLI version 1.0.0+.

For Python API documentation, see the main QCManyBody documentation.
