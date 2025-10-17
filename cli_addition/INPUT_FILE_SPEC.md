# QCManyBody CLI - Input File Specification

## Overview

The QCManyBody CLI accepts input files in YAML or JSON format. YAML is recommended for human readability, while JSON is useful for programmatic generation or integration with other tools.

## File Format

### Supported Formats
- **YAML** (`.yaml`, `.yml`) - Primary format, human-readable
- **JSON** (`.json`) - Alternative format, machine-readable

Format is auto-detected based on file extension and content.

## Schema Version

All input files should specify the schema version for compatibility:

```yaml
schema_version: "1.0"
```

## Top-Level Structure

```yaml
schema_name: "qcmanybody_cli_input"
schema_version: 1
molecule: {...}           # Required: Molecular system specification
calculation: {...}        # Required: Calculation parameters
bsse: {...}              # Optional: BSSE treatment (default: nocp)
manybody: {...}          # Optional: Many-body expansion parameters
output: {...}            # Optional: Output preferences
```

---

## Section: `molecule`

Specifies the molecular system to be calculated.

### Option 1: External File Reference

```yaml
molecule:
  source: xyz  # xyz, pdb, or qcschema
  file: path/to/molecule.xyz
  fragments: [[0,1,2], [3,4,5], [6,7,8]]  # Optional: Atom indices for each fragment
```

### Option 2: Inline Specification

```yaml
molecule:
  source: inline
  inline:
    symbols: [O, H, H, O, H, H, O, H, H]
    geometry: [
      [0.000, 0.000, 0.000],
      [0.758, 0.587, 0.000],
      [-0.758, 0.587, 0.000],
      [2.500, 0.000, 0.000],
      [3.258, 0.587, 0.000],
      [1.742, 0.587, 0.000],
      [5.000, 0.000, 0.000],
      [5.758, 0.587, 0.000],
      [4.242, 0.587, 0.000]
    ]
    fragments: [[0,1,2], [3,4,5], [6,7,8]]
    molecular_charge: 0.0
    molecular_multiplicity: 1
    units: angstrom  # or bohr
```

### Option 3: QCSchema Format

```yaml
molecule:
  source: inline
  format: qcschema
  symbols: [O, H, H, O, H, H, O, H, H]
  geometry: [  # In Bohr
    [0.000, 0.000, 0.000],
    [1.432, 1.109, 0.000],
    [-1.432, 1.109, 0.000],
    # ... more coordinates
  ]
  fragments: [[0,1,2], [3,4,5], [6,7,8]]
  molecular_charge: 0
  molecular_multiplicity: 1
```

### Fragment Specification

Fragments are specified as lists of atom indices (0-based):

```yaml
fragments: [[0,1,2], [3,4,5], [6,7,8]]  # 3 fragments
```

Or using fragment connectivity:

```yaml
fragments:
  - atoms: [0,1,2]
    charge: 0
    multiplicity: 1
  - atoms: [3,4,5]
    charge: 0
    multiplicity: 1
  - atoms: [6,7,8]
    charge: -1
    multiplicity: 1
```

---

## Section: `calculation`

Specifies the type of calculation and computational method.

### Single-Level Calculation

```yaml
calculation:
  driver: energy          # energy, gradient, hessian, properties
  method: mp2             # QC method (lowercase)
  basis: cc-pvdz          # Basis set (lowercase)
```

### Multi-Level Calculation

```yaml
calculation:
  driver: gradient
  levels:
    1: ccsd/cc-pvtz       # n-body level: method/basis
    2: mp2/cc-pvdz
    3: mp2/cc-pvdz
    # Higher levels use last specified level by default
```

**Notes:**
- `driver` must be one of: `energy`, `gradient`, `hessian`, `properties`
- Method and basis names are case-insensitive but stored as lowercase
- Multi-level syntax: `{n-body level}: {method}/{basis}`

---

## Section: `bsse`

Specifies Basis Set Superposition Error (BSSE) treatment.

### Simple BSSE Type

```yaml
bsse:
  type: cp  # nocp, cp, vmfc
```

### Multiple BSSE Types

```yaml
bsse:
  type: [cp, vmfc]  # Compute both CP and VMFC corrections
```

**Options:**
- `nocp`: No counterpoise correction
- `cp`: Standard counterpoise correction
- `vmfc`: Valiron-Mayer function counterpoise

---

## Section: `manybody`

Controls many-body expansion parameters.

```yaml
manybody:
  max_nbody: 3                    # Maximum n-body level (default: nfragments)
  return_total_data: true         # Return total energies (default: true)
  supersystem_ie_only: false      # Only compute supersystem IE (default: false)
  embedding_charges:              # Optional point charges
    - [1.0, 0.0, 0.0, 0.5]       # [x, y, z, charge] in Bohr
    - [2.0, 0.0, 0.0, -0.5]
```

**Parameters:**
- `max_nbody`: Maximum n-body level to compute (default: number of fragments)
- `return_total_data`: If true, return total energies; if false, only interaction energies
- `supersystem_ie_only`: If true, only compute supersystem interaction energy
- `embedding_charges`: List of point charges as [x, y, z, charge] in Bohr

---

## Section: `program`

Specifies the QC program and program-specific settings.

```yaml
program:
  name: psi4              # psi4, nwchem, cfour
  keywords:               # Program-specific keywords
    reference: rhf
    scf_type: df
    mp2_type: df
    freeze_core: true
  memory: 4GB            # Memory per process
```

**Supported Programs:**
- `psi4`: Psi4 quantum chemistry program
- `nwchem`: NWChem quantum chemistry program
- `cfour`: CFOUR quantum chemistry program

**Notes:**
- `keywords` are passed directly to the QC program
- Keyword names and values are program-specific
- Memory format: `{number}{unit}` where unit is B, KB, MB, GB, TB

---

**Note**: Parallel execution, resource management, and checkpoint/resume features are being developed in a separate branch and are not part of this CLI implementation.

---

## Section: `output`

Controls output format and verbosity.

```yaml
output:
  file: results.json             # Output file path (default: stdout)
  format: json                   # json, yaml, text
  verbosity: normal              # quiet, normal, verbose, debug
  save_intermediates: false      # Save intermediate results
  intermediate_dir: ./intermediates
```

**Format Options:**
- **`json`**: Machine-readable JSON (default)
- **`yaml`**: Human-readable YAML
- **`text`**: Human-readable text summary

**Verbosity Levels:**
- **`quiet`**: Only errors and final results
- **`normal`**: Standard output (default)
- **`verbose`**: Detailed progress information
- **`debug`**: Maximum detail for debugging

---

## Complete Examples

### Example 1: Basic Energy Calculation

```yaml
schema_name: "qcmanybody_cli_input"
schema_version: 1

molecule:
  source: xyz
  file: water_dimer.xyz
  fragments: [[0,1,2], [3,4,5]]

calculation:
  type: single
  single:
    driver: energy
    method: mp2
    basis: cc-pvdz
    program: psi4

bsse:
  type: [cp]

manybody:
  max_nbody: 2

output:
  file: results.json
  format: json
```

### Example 2: Multi-Level Gradient Calculation

```yaml
schema_name: "qcmanybody_cli_input"
schema_version: 1

molecule:
  source: xyz
  file: water_trimer.xyz
  fragments: [[0,1,2], [3,4,5], [6,7,8]]

calculation:
  type: multi
  multi:
    driver: gradient
    levels:
      1: {method: ccsd, basis: cc-pvtz, program: psi4}
      2: {method: mp2, basis: cc-pvdz, program: psi4}
      3: {method: mp2, basis: cc-pvdz, program: psi4}

bsse:
  type: [cp, vmfc]

manybody:
  max_nbody: 3
  return_total_data: true

output:
  file: results.yaml
  format: yaml
```

### Example 3: Inline Molecule with Embedding Charges

```yaml
schema_version: "1.0"

molecule:
  source: inline
  format: qcschema
  symbols: [Ne, Ne, Ne]
  geometry: [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 3.78],
    [0.0, 0.0, 7.56]
  ]
  fragments: [[0], [1], [2]]

calculation:
  driver: energy
  method: ccsd
  basis: cc-pvdz

bsse:
  type: cp

manybody:
  max_nbody: 3
  embedding_charges:
    - [5.0, 0.0, 0.0, 0.5]
    - [-5.0, 0.0, 0.0, -0.5]

program:
  name: nwchem
  keywords:
    ccsd:
      thresh: 1.0e-7

output:
  format: text
```

### Example 4: Multi-Level Energy Calculation

```yaml
schema_name: "qcmanybody_cli_input"
schema_version: 1

molecule:
  source: xyz
  file: large_cluster.xyz
  fragments: [[0,1,2,3,4], [5,6,7,8,9], [10,11,12,13,14], [15,16,17,18,19]]

calculation:
  type: multi
  multi:
    driver: energy
    levels:
      1: {method: ccsd, basis: aug-cc-pvtz, program: psi4}
      2: {method: mp2, basis: aug-cc-pvdz, program: psi4}
      3: {method: mp2, basis: cc-pvdz, program: psi4}
      4: {method: hf, basis: cc-pvdz, program: psi4}

bsse:
  type: [cp]

manybody:
  max_nbody: 4

output:
  file: results.json
  format: json
```

**Note:** Parallel execution and checkpointing features are being developed in a separate branch and are not yet available in this CLI version.

---

## Validation

The CLI will validate input files and provide clear error messages for:

1. **Schema errors**: Missing required fields, invalid field types
2. **Molecular errors**: Invalid geometry, fragment specification
3. **Method/basis errors**: Unsupported method or basis set
4. **Compatibility errors**: Incompatible option combinations
5. **Resource errors**: Invalid memory specifications, missing files

### Example Error Messages

```
Error in input.yaml:
  Section: molecule
  Issue: Fragment specification is incomplete
    Found 9 atoms but fragments only cover atoms [0, 1, 2, 3, 4, 5]
    Missing atoms: [6, 7, 8]
  Suggestion: Add fragment for atoms [6, 7, 8] or adjust atom count
```

```
Error in input.yaml:
  Section: calculation.levels
  Issue: Multi-level calculation requires at least 2 levels
    Found only 1 level defined
  Suggestion: Either define multiple levels or use single-level format:
    calculation:
      method: mp2
      basis: cc-pvdz
```

---

## Schema Validation

Input files are validated against a Pydantic schema. To see the complete schema:

```bash
qcmanybody validate --show-schema
```

To validate an input file without running:

```bash
qcmanybody validate input.yaml
```

---

## Conversion Between Formats

Convert between YAML and JSON:

```bash
# YAML to JSON
qcmanybody convert input.yaml output.json

# JSON to YAML
qcmanybody convert input.json output.yaml
```

---

## Environment Variable Substitution

Input files support environment variable substitution using `${VAR}` or `$VAR` syntax:

```yaml
molecule:
  source: xyz
  file: ${MOLECULE_DIR}/water_trimer.xyz

output:
  file: ${OUTPUT_DIR}/results_${RUN_ID}.json
```

---

## Comments

YAML files support comments using `#`:

```yaml
# This is a water trimer calculation
molecule:
  source: xyz
  file: water_trimer.xyz  # Path relative to input file
  fragments: [[0,1,2], [3,4,5], [6,7,8]]  # 3 water molecules
```

JSON does not support comments natively.

---

## Best Practices

1. **Use YAML for human-edited files** - More readable, supports comments
2. **Use JSON for programmatic generation** - Easier to generate from code
3. **Include schema_version** - Ensures compatibility with future versions
4. **Add comments in YAML** - Document your calculation setup
5. **Use relative paths** - Makes input files portable
6. **Start simple** - Begin with minimal required fields, add options as needed
7. **Validate before running** - Use `qcmanybody validate` to catch errors early
8. **Use dry-run mode** - Check execution plan before running expensive calculations

---

## Migration from Python Scripts

To convert existing Python scripts to input files:

1. Extract molecule specification
2. Extract calculation parameters
3. Map to input file schema
4. Test with `qcmanybody validate`
5. Run with `qcmanybody run`

Automated conversion tool:

```bash
qcmanybody convert script.py output.yaml --from python --to yaml
```

---

## Schema Versioning

The schema version follows semantic versioning:

- **Major version**: Breaking changes to schema
- **Minor version**: Backward-compatible additions
- **Patch version**: Bug fixes, clarifications

Current version: `1.0`

The CLI will warn if input file uses a different schema version and attempt to migrate if possible.
