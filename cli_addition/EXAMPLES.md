# QCManyBody CLI - Usage Examples

## Quick Start

### Installation

```bash
# Install QCManyBody with CLI support
pip install qcmanybody[cli]

# Or in development mode
cd QCManyBody
pip install -e ".[cli]"
```

### Verify Installation

```bash
# Check version
qcmanybody --version

# Show help
qcmanybody --help

# Show help for run command
qcmanybody run --help
```

---

## Basic Usage Examples

### Example 1: Simple Energy Calculation

**Input file** (`water_dimer_energy.yaml`):
```yaml
schema_version: "1.0"

molecule:
  source: file
  path: water_dimer.xyz
  fragments: [[0,1,2], [3,4,5]]

calculation:
  driver: energy
  method: mp2
  basis: cc-pvdz

bsse:
  type: cp

program:
  name: psi4
```

**Command:**
```bash
qcmanybody run water_dimer_energy.yaml -o results.json
```

**Output:** `results.json` contains the calculation results in JSON format.

---

### Example 2: Gradient Calculation

**Input file** (`water_dimer_gradient.yaml`):
```yaml
schema_version: "1.0"

molecule:
  source: file
  path: water_dimer.xyz
  fragments: [[0,1,2], [3,4,5]]

calculation:
  driver: gradient
  method: mp2
  basis: cc-pvdz

bsse:
  type: cp

program:
  name: psi4

output:
  format: yaml
```

**Command:**
```bash
qcmanybody run water_dimer_gradient.yaml -o results.yaml
```

---

### Example 3: Multi-Level Calculation

**Input file** (`water_trimer_multilevel.yaml`):
```yaml
schema_version: "1.0"

molecule:
  source: file
  path: water_trimer.xyz
  fragments: [[0,1,2], [3,4,5], [6,7,8]]

calculation:
  driver: gradient
  levels:
    1: ccsd/cc-pvtz
    2: mp2/cc-pvdz
    3: mp2/cc-pvdz

bsse:
  type: [cp, vmfc]

manybody:
  max_nbody: 3
  return_total_data: true

program:
  name: psi4

output:
  format: yaml
  verbosity: verbose
```

**Command:**
```bash
qcmanybody run water_trimer_multilevel.yaml -o results.yaml --verbose
```

---

## Planning and Validation

### Example 4: Show Execution Plan

Before running an expensive calculation, preview what will be computed:

```bash
qcmanybody plan water_trimer_multilevel.yaml --show-tasks
```

**Output:**
```
Execution Plan for: water_trimer_multilevel.yaml
================================================

Molecule: 3 fragments (9 atoms total)
Driver: gradient
Max n-body: 3
BSSE types: [cp, vmfc]

Tasks to be computed:
---------------------

1-body calculations (3 tasks):
  [1] CCSD/cc-pvtz: Fragment 1 in [1, 2, 3] (CP)
  [2] CCSD/cc-pvtz: Fragment 2 in [1, 2, 3] (CP)
  [3] CCSD/cc-pvtz: Fragment 3 in [1, 2, 3] (CP)
  ... (VMFC tasks)

2-body calculations (9 tasks):
  [4] MP2/cc-pvdz: Fragments [1, 2] in [1, 2, 3] (CP)
  [5] MP2/cc-pvdz: Fragments [1, 3] in [1, 2, 3] (CP)
  [6] MP2/cc-pvdz: Fragments [2, 3] in [1, 2, 3] (CP)
  ... (more CP and VMFC tasks)

3-body calculations (2 tasks):
  [10] MP2/cc-pvdz: Fragments [1, 2, 3] in [1, 2, 3] (CP)
  [11] MP2/cc-pvdz: Fragments [1, 2, 3] in [1, 2, 3] (VMFC)

Total tasks: 24
Estimated time: ~45 minutes (assuming 2 min per task)
Estimated memory: ~4 GB per task
```

---

### Example 5: Validate Input File

Check if input file is valid before running:

```bash
qcmanybody validate water_trimer_multilevel.yaml
```

**Output (success):**
```
✓ Input file is valid
  Schema version: 1.0
  Molecule: 9 atoms, 3 fragments
  Calculation: gradient, multi-level (3 levels)
  Program: psi4
```

**Output (with errors):**
```
✗ Input file has errors:

Error in molecule section:
  Line 5: Fragment specification incomplete
    Found 9 atoms but fragments only cover 6 atoms
    Missing atoms: [6, 7, 8]
  Suggestion: Add fragment [[6, 7, 8]]

Error in calculation section:
  Line 12: Multi-level requires at least 2 levels
    Found only 1 level defined
```

---

## Advanced Features

**Note**: Parallel execution, checkpointing, and resource management features are being developed in a separate branch and are not part of this CLI implementation.

---

### Example 6: Verbose Logging

Run with detailed logging:

```bash
qcmanybody run water_trimer.yaml \
  --verbose \
  -o results.json
```

This will print detailed progress information during the calculation.

---

## Working with Different Molecule Formats

### Example 7: XYZ File

**water_dimer.xyz:**
```
6
Water dimer
O   0.000   0.000   0.000
H   0.758   0.587   0.000
H  -0.758   0.587   0.000
O   2.500   0.000   0.000
H   3.258   0.587   0.000
H   1.742   0.587   0.000
```

**Input file:**
```yaml
molecule:
  source: file
  path: water_dimer.xyz
  format: xyz
  fragments: [[0,1,2], [3,4,5]]
```

---

### Example 8: Inline Molecule Definition

```yaml
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
  molecular_charge: 0
  molecular_multiplicity: 1
```

---

### Example 9: PDB File (Future)

```yaml
molecule:
  source: file
  path: protein_fragment.pdb
  format: pdb
  fragments: [[0,1,2,...,50], [51,52,...,100], [101,102,...,150]]
```

---

## Different BSSE Treatments

### Example 10: No Counterpoise

```yaml
bsse:
  type: nocp

calculation:
  driver: energy
  method: mp2
  basis: cc-pvdz
```

---

### Example 11: Multiple BSSE Methods

Compute both CP and VMFC corrections:

```yaml
bsse:
  type: [cp, vmfc]

calculation:
  driver: gradient
  method: mp2
  basis: cc-pvdz
```

The output will contain results for both methods.

---

## Output Formats

### Example 12: JSON Output (Default)

```bash
qcmanybody run input.yaml -o results.json
```

**results.json:**
```json
{
  "schema_version": "1.0",
  "molecule": {...},
  "calculation": {...},
  "results": {
    "cp": {
      "energy": -152.123456,
      "gradient": [[...], [...], ...],
      "1b_energy": {...},
      "2b_energy": {...}
    }
  },
  "timing": {
    "total": 123.45,
    "tasks": [...]
  }
}
```

---

### Example 13: YAML Output

```bash
qcmanybody run input.yaml -o results.yaml --format yaml
```

**results.yaml:**
```yaml
schema_version: "1.0"
molecule:
  # ...
calculation:
  # ...
results:
  cp:
    energy: -152.123456
    gradient:
      - [0.0, 0.0, 0.001]
      - [0.0, 0.0, -0.001]
      # ...
```

---

### Example 14: Text Summary

```bash
qcmanybody run input.yaml --format text
```

**Output:**
```
QCManyBody Calculation Results
==============================

Molecule: 6 atoms, 2 fragments
Method: MP2/cc-pVDZ
BSSE: CP
Driver: gradient

Results (CP):
-------------
Total Energy:    -152.123456 Eh
Interaction Energy: -0.005432 Eh (-3.41 kcal/mol)

1-body Energy:   -152.118024 Eh
2-body Energy:   -0.005432 Eh

Gradient (Eh/Bohr):
Fragment 1:
  O:  [  0.000000,  0.000000,  0.001234 ]
  H:  [  0.000567, -0.000123, -0.000617 ]
  H:  [ -0.000567, -0.000123, -0.000617 ]
Fragment 2:
  O:  [  0.000000,  0.000000, -0.001234 ]
  H:  [ -0.000567,  0.000123,  0.000617 ]
  H:  [  0.000567,  0.000123,  0.000617 ]

Computation Time: 123.45 seconds
Tasks Completed: 12/12
```

---

## Conversion Between Formats

### Example 15: YAML to JSON

```bash
qcmanybody convert input.yaml output.json
```

---

### Example 16: JSON to YAML

```bash
qcmanybody convert input.json output.yaml
```

---

## Integration with Shell Scripts

### Example 17: Batch Processing

**run_batch.sh:**
```bash
#!/bin/bash

# Run multiple calculations
for mol in molecules/*.xyz; do
  name=$(basename "$mol" .xyz)

  # Create input file
  cat > "${name}_input.yaml" <<EOF
schema_version: "1.0"
molecule:
  source: file
  path: $mol
  fragments: [[0,1,2], [3,4,5]]
calculation:
  driver: energy
  method: mp2
  basis: cc-pvdz
bsse:
  type: cp
program:
  name: psi4
output:
  path: ${name}_results.json
EOF

  # Run calculation
  qcmanybody run "${name}_input.yaml" --verbose
done
```

---

### Example 18: Environment Variables

**Input file** (`template.yaml`):
```yaml
schema_version: "1.0"

molecule:
  source: file
  path: ${MOLECULE_PATH}

calculation:
  driver: energy
  method: ${METHOD}
  basis: ${BASIS}

program:
  name: ${PROGRAM}

output:
  path: ${OUTPUT_PATH}
```

**Run script:**
```bash
#!/bin/bash

export MOLECULE_PATH="water_dimer.xyz"
export METHOD="mp2"
export BASIS="cc-pvdz"
export PROGRAM="psi4"
export OUTPUT_PATH="results.json"

qcmanybody run template.yaml
```

---

## Tips and Best Practices

1. **Always validate before running:**
   ```bash
   qcmanybody validate input.yaml && qcmanybody run input.yaml
   ```

2. **Use plan mode for expensive calculations:**
   ```bash
   qcmanybody plan input.yaml --show-tasks
   ```

3. **Use verbose mode for debugging:**
   ```bash
   qcmanybody run input.yaml --verbose
   ```

4. **Redirect output for batch jobs:**
   ```bash
   qcmanybody run input.yaml > stdout.log 2> stderr.log
   ```

5. **Monitor progress in real-time:**
   ```bash
   qcmanybody run input.yaml --verbose 2>&1 | tee progress.log
   ```
