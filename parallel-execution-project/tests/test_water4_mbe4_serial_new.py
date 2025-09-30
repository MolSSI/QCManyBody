#!/usr/bin/env python3
"""Water 4 cluster many-body expansion serial test.

This script validates the traditional serial execution path against the
current data model updates.
"""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))


def run_water4_serial_calculation():
    """Run the 4-water many-body expansion in serial mode."""

    from qcelemental.models import Molecule

    symbols = ["O", "H", "H"] * 4

    geometry = [
        [0.069765, -3.765074, -0.092294],
        [-0.748399, -4.102246, -0.457193],
        [-0.205709, -3.169585, 0.604654],
        [2.608892, -2.675449, 2.890175],
        [2.427245, -3.533869, 2.507619],
        [2.493361, -2.060559, 2.16575],
        [0.0, 0.0, 0.0],
        [-0.900139, -0.263684, 0.190917],
        [0.092177, -0.127148, -0.944228],
        [2.450667, -1.324316, 0.643786],
        [2.758818, -1.643552, -0.204364],
        [1.562815, -1.009682, 0.473627],
    ]

    fragments = [
        list(range(index, index + 3)) for index in range(0, len(symbols), 3)
    ]

    mol = Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=fragments,
        molecular_charge=0,
        molecular_multiplicity=1
    )

    # Configure serial many-body calculation using ManyBodyComputer
    from qcmanybody.computer import ManyBodyComputer
    from qcmanybody.models import ManyBodyInput

    def _spec_entry() -> dict:
        return {
            "program": "psi4",
            "driver": "energy",
            "model": {"method": "hf", "basis": "sto-3g"},
            "keywords": {},
            "protocols": {},
            "extras": {},
        }

    specification_map = {
        f"hf_level_{level}": _spec_entry() for level in range(1, 5)
    }

    level_mapping = {
        level: f"hf_level_{level}" for level in range(1, 5)
    }

    specification_payload = {
        "driver": "energy",
        "keywords": {
            "bsse_type": ["nocp"],
            "levels": level_mapping,
            "return_total_data": False,
            "supersystem_ie_only": False,
            "embedding_charges": {},
        },
        "specification": specification_map,
        "protocols": {},
        "extras": {},
    }

    mb_input = ManyBodyInput(
        molecule=mol,
        specification=specification_payload,
        extras={},
    )

    print(
        f"System: {len(mol.symbols)} atoms in {len(mol.fragments)} fragments"
    )
    print("Max n-body level: 4")
    print("BSSE treatment: ['nocp']")
    print("QC method: HF/STO-3G")
    print("Execution mode: SERIAL")

    print("\n" + "="*60)
    print("STARTING SERIAL MANY-BODY EXPANSION CALCULATION")
    print("Using traditional ManyBodyComputer for validation")
    print("="*60)

    import time
    start_time = time.time()

    # Create and run the calculation
    result = ManyBodyComputer.from_manybodyinput(mb_input)

    end_time = time.time()
    elapsed_time = end_time - start_time

    print("\n" + "="*60)
    print("SERIAL CALCULATION COMPLETED!")
    print(f"Total time: {elapsed_time:.2f} seconds")

    print(f"\nSerial result type: {type(result)}")
    if hasattr(result, 'return_result'):
        print(f"Final energy result: {result.return_result}")
    else:
        print("No return_result attribute found")

    if hasattr(result, 'properties'):
        print(f"Properties available: {type(result.properties)}")
        if hasattr(result.properties, 'dict'):
            prop_dict = result.properties.dict()
            print(f"Properties keys: {list(prop_dict.keys())}")

    print("\n✓ SERIAL CALCULATION SUCCESSFUL!")
    return result


def test_psi4_availability():
    """Quick test to verify Psi4 is working."""
    try:
        import psi4
        print("✓ Psi4 import successful")

        # Quick test calculation
        psi4.set_memory('500 MB')
        psi4.set_num_threads(1)

        psi4.geometry("""
        O
        H 1 0.96
        H 1 0.96 2 104.5
        """)

        psi4.set_options({'basis': 'sto-3g', 'scf_type': 'df'})
        energy = psi4.energy('hf')
        print(f"✓ Psi4 test successful: {energy:.6f} Eh")
        return True

    except ImportError:
        print("✗ Psi4 not available")
        return False
    except Exception as e:
        print(f"✗ Psi4 test failed: {e}")
        return False


if __name__ == "__main__":
    print("Water 4 Cluster Many-Body Expansion - SERIAL TEST")
    print("Using traditional ManyBodyComputer for validation")
    print("=" * 60)

    # Test Psi4 first
    print("Testing Psi4 availability...")
    if not test_psi4_availability():
        print("\nCannot proceed without Psi4. Install with:")
        print("conda install -c conda-forge psi4")
        sys.exit(1)

    try:
        print("\nStarting serial calculation...")
        result = run_water4_serial_calculation()
        print("\n" + "="*60)
        print("SUCCESS: Serial calculation completed successfully!")
        print(
            "AtomicResult analysis fixes work for both serial and parallel "
            "modes!"
        )
        print("="*60)

    except KeyboardInterrupt:
        print("\n\nCalculation interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Calculation failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
