# Units Conversion Bug Fix Summary

## Date: January 5, 2026

## Critical Bug Identified

A systematic units conversion bug was found in QCManyBody manuscript calculation scripts, identical to the bug previously fixed in the HMBE package.

### Root Cause

**QCElemental Molecule Default Behavior:**
- `Molecule()` constructor expects coordinates in **Bohr** (atomic units) by default
- `Molecule.from_data()` with `"units angstrom"` properly converts to Bohr internally
- XYZ files use **Angstrom** by convention

**The Bug:**
Scripts were reading XYZ coordinates (Angstrom) and passing them directly to `Molecule()` constructor, which interpreted them as Bohr. This caused all molecular geometries to be **shrunk by ~47%** (the Bohr→Angstrom conversion factor of 0.529177).

### Impact

**Before Fix:**
- O-H bond length: 0.51 Å (wrong - shrunk by 47%)
- Total energies significantly incorrect
- Many-body contributions completely wrong

**After Fix:**
- O-H bond length: 0.97 Å (correct - expected ~0.96 Å)
- Geometries match physical reality
- Energies and MBE contributions accurate

## Files Fixed

### 1. water16.*_qcmanybody.py Scripts (22 files)

**Location:** `parallel-execution-project/manuscript_calculations/ref_conformations/`

**Issue:** `load_xyz()` function read XYZ files (Angstrom) but created Molecule objects without unit conversion.

**Fix:** Replaced `load_xyz()` and `build_fragments()` with new `load_xyz_to_molecule()` function that:
1. Reads XYZ file coordinates
2. Uses `Molecule.from_data()` with explicit `"units angstrom"` declaration
3. Returns properly configured Molecule object

**Files affected:**
- `water16.000_qcmanybody.py` through `water16.010_qcmanybody.py` (11 files)
- `mbe4/water16.000_qcmanybody.py` through `mbe4/water16.010_qcmanybody.py` (11 files)

### 2. water_*_qcmanybody.py Scripts (14 files)

**Location:** `parallel-execution-project/manuscript_calculations/`

**Issue:** Hardcoded GEOMETRY in Angstrom passed directly to Molecule() constructor.

**Fix:** Modified `create_molecule()` to convert Angstrom → Bohr:
```python
from qcelemental import constants

geometry_bohr = [
    [x / constants.bohr2angstroms, y / constants.bohr2angstroms, z / constants.bohr2angstroms]
    for x, y, z in GEOMETRY
]

return Molecule(
    symbols=SYMBOLS,
    geometry=geometry_bohr,  # Now in Bohr
    fragments=FRAGMENTS,
    molecular_charge=TOTAL_CHARGE,
    molecular_multiplicity=TOTAL_MULTIPLICITY,
)
```

**Files affected:**
- `water_000_qcmanybody.py`
- `water_000_qcmanybody_hf.py`
- `water_000_qcmanybody_ccsdt.py`
- `water_001_qcmanybody.py` through `water_010_qcmanybody.py`
- `water_ref1_qcmanybody.py`

### 3. Test Files (1 file)

**Location:** `parallel-execution-project/tests/`

**Issue:** `test_water4_mbe4_multiprocessing.py` had hardcoded geometry in Angstrom.

**Fix:** Added Angstrom → Bohr conversion in `build_water4_molecule()` function.

**File affected:**
- `test_water4_mbe4_multiprocessing.py`

### 4. Template Generator (1 file)

**Location:** `parallel-execution-project/manuscript_calculations/`

**Issue:** `hmbe_to_qcmanybody.py` template generated scripts with the bug.

**Fix:** Updated template's `create_molecule()` function to include Angstrom → Bohr conversion with documentation explaining why.

**File affected:**
- `hmbe_to_qcmanybody.py`

## QCManyBody Core Library Status

✅ **NO BUGS FOUND** in the core QCManyBody library (`qcmanybody/` directory)

The core library correctly handles units:
- All test molecules use `Molecule.from_data()` with explicit `"units bohr"`
- Direct `Molecule()` calls in tests use Bohr values
- Fragment extraction via `get_fragment()` preserves units correctly
- No incorrect Bohr↔Angstrom conversions in core code

## Verification

Fixed code was validated with test calculations:

```python
# XYZ file coordinates (Angstrom): [0.691332, -2.695712, 1.453397]
# Expected O-H distance: ~0.96 Å

# BEFORE FIX (treating Angstrom as Bohr):
O-H distance: 0.5127 Å  ✗ WRONG

# AFTER FIX (converting Angstrom to Bohr):
O-H distance: 0.9689 Å  ✓ CORRECT
```

## Best Practices Going Forward

1. **Always specify units explicitly** when creating Molecule objects
2. **Use `Molecule.from_data()` with `"units angstrom"`** when reading from XYZ files
3. **Document units** at every stage: parser output, storage format, QC backend input
4. **Add unit tests** to verify geometry handling from various input formats
5. **When using `Molecule()` constructor directly**, remember it expects **Bohr** by default

## References

- Similar bug was previously fixed in HMBE package (December 2025)
- QCElemental documentation: https://github.com/MolSSI/QCElemental
- Bohr to Angstrom conversion factor: 0.52917721067 (from `qcelemental.constants.bohr2angstroms`)
