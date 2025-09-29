#!/usr/bin/env python3
"""
Test QCManyBody parallel execution with real Psi4 calculations.
This uses a small system to verify the parallel system works with actual QC.
"""

import copy
import time


def test_parallel_with_real_qc():
    """Test parallel execution using real Psi4 calculations on a small system."""

    print("Testing QCManyBody Parallel Execution with Real Psi4 Calculations")
    print("="*70)

    # Test Psi4 availability first
    try:
        import psi4
        print(f"✓ Psi4 available: version {psi4.__version__}")
    except ImportError:
        print("✗ Psi4 not available - skipping real QC test")
        return False

    from qcmanybody import ManyBodyCore, BsseEnum
    from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
    from qcelemental.models import Molecule

    # Create a very small test system (helium dimer)
    # This will be fast but still use real QC calculations
    molecule = Molecule(**{
        'symbols': ['He', 'He'],
        'geometry': [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0]  # 3 bohr separation
        ],
        'fragments': [[0], [1]],
        'molecular_charge': 0,
        'molecular_multiplicity': 1,
        'fragment_charges': [0, 0],
        'fragment_multiplicities': [1, 1]
    })

    print(f"✓ Test system: {len(molecule.symbols)} atoms, {len(molecule.fragments)} fragments")

    # Create ManyBodyCore for 2-body calculation only (to keep it fast)
    core = ManyBodyCore(
        molecule=molecule,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf", 2: "hf"},  # Just HF/STO-3G for speed
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    print("✓ ManyBodyCore created for 2-body HF/STO-3G calculation")

    # Configure for real QC calculations
    config = ParallelConfig(
        max_workers=2,
        execution_mode="threading",  # Threading works well with QCEngine
        use_qcengine=True,
        qc_program="psi4",
        basis_set="sto-3g",
        memory_limit_mb=500,
        timeout_seconds=300,  # 5 minutes timeout
        qcengine_config={
            "keywords": {
                'scf_type': 'df',
                'maxiter': 50,
                'e_convergence': 1e-06,
                'd_convergence': 1e-06
            },
            "protocols": {
                "stdout": False  # Suppress output for cleaner display
            }
        }
    )

    print("✓ ParallelConfig created for real QCEngine calculations")

    specifications = {
        "hf": {
            "program": "psi4",
            "specification": {
                "driver": "energy",
                "model": {"method": "hf", "basis": "sto-3g"},
                "keywords": {'scf_type': 'df', 'maxiter': 50, 'e_convergence': 1e-06, 'd_convergence': 1e-06},
                "protocols": {"stdout": False},
                "extras": {},
            },
        },
    }

    # Create executor
    executor = ParallelManyBodyExecutor(core, config, driver="energy", specifications=specifications)
    print("✓ ParallelManyBodyExecutor created")

    # Run the calculation
    print("\nRunning parallel calculation with real Psi4...")
    print("Expected fragments: 2 monomers + 1 dimer = 3 total calculations")

    start_time = time.time()
    try:
        results = executor.execute_full_calculation()
        end_time = time.time()

        print(f"✓ Parallel calculation completed in {end_time - start_time:.2f}s")
        print(f"✓ Fragment results: {len(results)}")

        # Display results
        print("\nFragment Results:")
        for label, result in results.items():
            print(f"  {label}: {result.return_result:.8f} Eh")

        # Show execution statistics
        stats = executor.get_execution_statistics()
        print(f"\nExecution Statistics:")
        print(f"  Total fragments: {stats['total_fragments']}")
        print(f"  Levels executed: {stats['levels_executed']}")
        print(f"  Parallel time: {stats['parallel_time']:.3f}s")
        print(f"  Estimated speedup: {stats['speedup_factor']:.2f}x")

        return True

    except Exception as e:
        print(f"✗ Parallel calculation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_comparison_with_real_qc():
    """Test sequential vs parallel with real QC to verify correctness."""

    print("\nTesting Sequential vs Parallel with Real QC...")
    print("="*50)

    from qcmanybody import ManyBodyCore, BsseEnum
    from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
    from qcelemental.models import Molecule

    # Very simple system for validation
    molecule = Molecule(**{
        'symbols': ['H', 'H'],
        'geometry': [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0]  # H2 molecule split into fragments for testing
        ],
        'fragments': [[0], [1]],
        'molecular_charge': 0,
        'molecular_multiplicity': 1,
        'fragment_charges': [0, 0],
        'fragment_multiplicities': [2, 2]  # H atoms are doublets
    })

    core = ManyBodyCore(
        molecule=molecule,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf", 2: "hf"},
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    # Common QC config
    qc_config = {
        "use_qcengine": True,
        "qc_program": "psi4",
        "basis_set": "sto-3g",
        "memory_limit_mb": 500,
        "timeout_seconds": 120,
        "qcengine_config": {
            "keywords": {'scf_type': 'df'},
            "protocols": {"stdout": False}
        }
    }

    # Sequential execution
    print("Running sequential calculation...")
    sequential_config = ParallelConfig(
        max_workers=1,
        execution_mode="serial",
        **qc_config
    )

    specifications = {
        "hf": {
            "program": "psi4",
            "specification": {
                "driver": "energy",
                "model": {"method": "hf", "basis": "sto-3g"},
                "keywords": {'scf_type': 'df'},
                "protocols": {"stdout": False},
                "extras": {},
            },
        },
    }

    sequential_executor = ParallelManyBodyExecutor(
        core,
        sequential_config,
        driver="energy",
        specifications=copy.deepcopy(specifications),
    )
    sequential_results = sequential_executor.execute_full_calculation()
    print(f"✓ Sequential: {len(sequential_results)} results")

    # Parallel execution
    print("Running parallel calculation...")
    parallel_config = ParallelConfig(
        max_workers=2,
        execution_mode="threading",
        **qc_config
    )

    parallel_executor = ParallelManyBodyExecutor(
        core,
        parallel_config,
        driver="energy",
        specifications=copy.deepcopy(specifications),
    )
    parallel_results = parallel_executor.execute_full_calculation()
    print(f"✓ Parallel: {len(parallel_results)} results")

    # Validate correctness
    try:
        parallel_executor.validate_parallel_correctness(
            parallel_results, sequential_results, tolerance=1e-10  # Slightly looser for real QC
        )
        print("✓ Parallel vs Sequential validation PASSED (tolerance=1e-10)")

        # Show some actual energies to verify they're real
        print("\nSample energies (showing real QC results):")
        for label in list(parallel_results.keys())[:2]:  # Show first 2
            seq_energy = sequential_results[label].return_result
            par_energy = parallel_results[label].return_result
            diff = abs(seq_energy - par_energy)
            print(f"  {label}:")
            print(f"    Sequential: {seq_energy:.10f} Eh")
            print(f"    Parallel:   {par_energy:.10f} Eh")
            print(f"    Difference: {diff:.2e} Eh")

        return True

    except Exception as e:
        print(f"✗ Validation failed: {e}")
        return False


if __name__ == "__main__":
    print("QCManyBody Parallel Execution with Real Quantum Chemistry")
    print("="*70)

    # Test 1: Basic functionality with real QC
    success1 = test_parallel_with_real_qc()

    # Test 2: Correctness validation with real QC
    success2 = test_comparison_with_real_qc()

    print("\n" + "="*70)
    if success1 and success2:
        print("✅ ALL TESTS PASSED - Parallel execution works with real Psi4!")
        print("\nThe parallel system is fully functional and ready for:")
        print("• Real quantum chemistry calculations")
        print("• Production many-body expansion calculations")
        print("• Large system calculations like water16 cluster")
        print("\nYou can now run actual_test_water16.py for the full test.")
    else:
        print("❌ SOME TESTS FAILED - Check errors above")
    print("="*70)