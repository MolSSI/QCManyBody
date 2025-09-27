#!/usr/bin/env python3
"""
Quick validation test for multiprocessing serialization fixes.

This test validates that P1A-002 serialization fixes work by testing
fragment execution with multiprocessing mode, without hitting the
analysis phase bug (documented as Issue #001).
"""

import sys
import os
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)

def test_multiprocessing_fragment_execution():
    """Test that fragment execution works with multiprocessing (no analysis phase)."""

    from qcelemental.models import Molecule
    from qcmanybody import ManyBodyCore, BsseEnum
    from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig

    # Simple water dimer for quick testing
    mol = Molecule.from_data("""
    O 0.0000000000  0.0000000000 -0.0657755706
    H 0.0000000000 -0.7590619907  0.5219530189
    H 0.0000000000  0.7590619907  0.5219530189
    O 2.5000000000  0.0000000000  0.0000000000
    H 3.2570000000  0.5860000000  0.0000000000
    H 1.7430000000  0.5860000000  0.0000000000
    """, molecular_charge=0, molecular_multiplicity=1, fragments=[[0, 1, 2], [3, 4, 5]])

    # Create ManyBodyCore with minimal setup
    core = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf", 2: "hf"},  # Only 1-body and 2-body
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    # Configure multiprocessing execution
    parallel_config = ParallelConfig(
        max_workers=2,
        execution_mode="multiprocessing",  # Testing multiprocessing serialization
        use_qcengine=True,
        qc_program="psi4",
        basis_set="sto-3g",  # Minimal basis for speed
        memory_limit_mb=500,
        timeout_seconds=60,
        qcengine_config={
            "keywords": {
                'scf_type': 'df',
                'maxiter': 50,
                'e_convergence': 1e-06,
                'd_convergence': 1e-06
            },
            "protocols": {
                "stdout": False
            }
        }
    )

    print("Testing multiprocessing fragment execution...")
    print(f"System: {len(mol.symbols)} atoms, {len(mol.fragments)} fragments")
    print(f"Mode: multiprocessing with {parallel_config.max_workers} workers")

    # Create executor and execute fragments (but skip analysis)
    executor = ParallelManyBodyExecutor(core, parallel_config)

    print("Executing fragments with multiprocessing...")
    fragment_results = executor.execute_full_calculation()

    # Get execution statistics without triggering analysis phase
    stats = executor.get_execution_statistics()

    print(f"\n✓ Multiprocessing execution successful!")
    print(f"  Total fragments executed: {stats['total_fragments']}")
    print(f"  N-body levels processed: {stats['levels_executed']}")
    print(f"  Parallel execution time: {stats['parallel_time']:.2f}s")
    print(f"  Fragment results type: {type(fragment_results)}")
    print(f"  Fragment results keys: {list(fragment_results.keys()) if hasattr(fragment_results, 'keys') else 'N/A'}")

    # Validate that we have the expected number of fragments
    expected_fragments = 3  # 2 monomers + 1 dimer
    if len(fragment_results) == expected_fragments:
        print(f"✓ Correct number of fragments: {expected_fragments}")
    else:
        print(f"✗ Fragment count mismatch: got {len(fragment_results)}, expected {expected_fragments}")
        return False

    print(f"\n✓ P1A-002 VALIDATION SUCCESSFUL: Multiprocessing serialization works!")
    return True

if __name__ == "__main__":
    try:
        success = test_multiprocessing_fragment_execution()
        if success:
            print("\n" + "="*60)
            print("SUCCESS: Multiprocessing serialization fixes validated!")
            print("P1A-002 implementation complete!")
            print("="*60)
        else:
            sys.exit(1)
    except Exception as e:
        print(f"\n✗ Validation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)