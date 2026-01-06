#!/usr/bin/env python3
"""
Simple test to verify QCEngine works properly with our setup.
"""

def test_qcengine_basic():
    """Test basic QCEngine functionality."""

    try:
        import qcengine as qcng
        from qcelemental.models import AtomicInput, Molecule

        print("Testing QCEngine basic functionality...")

        # Simple helium atom
        molecule = Molecule(**{
            'symbols': ['He'],
            'geometry': [[0.0, 0.0, 0.0]],
            'molecular_charge': 0,
            'molecular_multiplicity': 1
        })

        print("✓ Molecule created")

        # Create atomic input
        atomic_input = AtomicInput(
            molecule=molecule,
            driver="energy",
            model={"method": "hf", "basis": "sto-3g"},
            keywords={"scf_type": "df"}
        )

        print("✓ AtomicInput created")

        # Run calculation
        result = qcng.compute(
            atomic_input,
            "psi4",
            raise_error=True,
            task_config={"memory": 0.5, "ncores": 1}  # 0.5 GB
        )

        print(f"✓ QCEngine calculation successful: {result.return_result:.8f} Eh")
        print(f"✓ Success status: {result.success}")

        return True

    except Exception as e:
        print(f"✗ QCEngine test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    test_qcengine_basic()