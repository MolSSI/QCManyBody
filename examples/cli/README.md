# QCManyBody CLI Examples

This directory contains example input files for the QCManyBody command-line interface.

## Examples

### 01_basic_energy - Basic Energy Calculation
**Files**: `01_basic_energy.json`, `01_basic_energy.yaml`

The simplest possible many-body calculation: 3 helium atoms with MP2/cc-pVDZ.
Demonstrates:
- Inline molecule specification
- Single-level calculation
- Counterpoise correction
- JSON and YAML formats

**Usage**:
```bash
qcmanybody run 01_basic_energy.json
# or
qcmanybody run 01_basic_energy.yaml
```

### 02_gradient - Gradient Calculation
**File**: `02_gradient.json`

Calculates gradients for geometry optimization.
Demonstrates:
- Gradient driver
- Multiple BSSE types (CP and NoCP)
- Returning total data (not just interaction)

**Usage**:
```bash
qcmanybody run 02_gradient.json
```

### 03_multilevel - Multi-Level Calculation
**File**: `03_multilevel.json`

Uses different methods at different n-body levels (high accuracy for 1-body, lower for higher bodies).
Demonstrates:
- Multi-level specification
- Different programs for different levels
- Program-specific keywords
- VMFC correction

**Usage**:
```bash
qcmanybody run 03_multilevel.json
```

### 04_from_xyz - Load from XYZ File
**Files**: `04_from_xyz.json`, `water_dimer.xyz`

Loads molecular structure from an external XYZ file.
Demonstrates:
- External molecule file
- Fragment specification for file-based molecules
- Text output format
- Water dimer interaction

**Usage**:
```bash
cd examples/cli
qcmanybody run 04_from_xyz.json
```

## Input File Format

### JSON Format
- Requires no additional dependencies
- Strict syntax (no comments, trailing commas forbidden)
- Best for programmatic generation

Example:
```json
{
  "schema_name": "qcmanybody_cli_input",
  "schema_version": 1,
  "molecule": { ... },
  "calculation": { ... }
}
```

### YAML Format
- Requires PyYAML: `pip install qcmanybody[cli]`
- Supports comments
- More human-readable
- Best for manual editing

Example:
```yaml
schema_name: qcmanybody_cli_input
schema_version: 1

molecule:  # Molecule specification
  source: inline
  ...
```

## Running Examples

### Prerequisites
- QCManyBody installed: `pip install qcmanybody[standard]`
- QC program (Psi4, NWChem, or CFOUR) available
- For YAML examples: `pip install qcmanybody[cli]`

### Basic Usage
```bash
# Run a calculation
qcmanybody run input.json

# Validate input without running
qcmanybody validate input.json

# Show execution plan
qcmanybody plan input.json

# Convert between formats
qcmanybody convert input.json input.yaml
```

### With Options
```bash
# Save output to file
qcmanybody run input.json -o results.json

# Use multiple processors
qcmanybody run input.json --parallel multiprocessing --nproc 4

# Verbose output
qcmanybody -vv run input.json

# Quiet mode (errors only)
qcmanybody -q run input.json
```

## Creating Your Own Input Files

1. Start with the simplest example (01_basic_energy.json)
2. Modify the molecule, method, and options
3. Validate before running:
   ```bash
   qcmanybody validate my_input.json
   ```
4. Check the execution plan:
   ```bash
   qcmanybody plan my_input.json --show-tasks
   ```
5. Run the calculation:
   ```bash
   qcmanybody run my_input.json -o results.json
   ```

## Common Patterns

### Molecule Specification
```json
// Inline (for small molecules)
"molecule": {
  "source": "inline",
  "inline": {
    "symbols": ["He", "He"],
    "geometry": [[0,0,0], [0,0,2]],
    "fragments": [[0], [1]]
  }
}

// From XYZ file
"molecule": {
  "source": "xyz",
  "file": "molecule.xyz",
  "fragments": [[0,1,2], [3,4,5]]
}
```

### Single vs Multi-Level
```json
// Single level (same method for all n-body)
"calculation": {
  "type": "single",
  "single": {
    "driver": "energy",
    "method": "mp2",
    "basis": "cc-pvdz",
    "program": "psi4"
  }
}

// Multi-level (different methods)
"calculation": {
  "type": "multi",
  "multi": {
    "driver": "energy",
    "levels": {
      "1": {"method": "ccsd(t)", "basis": "cc-pvtz", "program": "cfour"},
      "2": {"method": "mp2", "basis": "cc-pvdz", "program": "psi4"}
    }
  }
}
```

### BSSE Corrections
```json
// Single BSSE type
"bsse": {
  "type": ["cp"]
}

// Multiple BSSE types (first is returned)
"bsse": {
  "type": ["cp", "nocp", "vmfc"]
}
```

## Troubleshooting

### Validation Errors
If validation fails, check:
- Required fields are present
- Field types are correct (numbers vs strings)
- Fragment indices are 0-based
- All atoms are assigned to fragments

### Runtime Errors
If calculation fails:
- Check QC program is installed and in PATH
- Verify method/basis compatibility with program
- Check molecule geometry is reasonable
- Try with simpler method first (e.g., SCF/STO-3G)

### Performance
For large calculations:
- Use `--parallel multiprocessing --nproc N`
- Set `max_nbody` to limit n-body levels
- Use `supersystem_ie_only: true` to skip intermediate bodies
- Consider multi-level: expensive method only for 1-body

## More Information

- Full documentation: https://molssi.github.io/QCManyBody/
- Input schema reference: `qcmanybody validate --show-schema`
- GitHub issues: https://github.com/MolSSI/QCManyBody/issues
