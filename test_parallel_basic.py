#!/usr/bin/env python3
"""
Quick test to verify the parallel execution system is working.
"""

import copy

def test_parallel_system():
    """Test basic parallel execution functionality."""

    # Import required modules
    from qcmanybody import ManyBodyCore, BsseEnum
    from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
    from qcelemental.models import Molecule

    print("Testing QCManyBody Parallel Execution System...")
    print("="*50)

    # Create a simple test molecule (water dimer)
    molecule = Molecule(**{
        'symbols': ['O', 'H', 'H', 'O', 'H', 'H'],
        'geometry': [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  1.43233673, -0.96104039],
            [ 0.0, -1.43233673, -0.96104039],
            [ 3.0,  0.0,  0.0],
            [ 3.0,  1.43233673, -0.96104039],
            [ 3.0, -1.43233673, -0.96104039],
        ],
        'fragments': [[0, 1, 2], [3, 4, 5]],
        'molecular_charge': 0.0,
        'molecular_multiplicity': 1
    })

    print(f"✓ Test molecule created: {len(molecule.symbols)} atoms, {len(molecule.fragments)} fragments")

    # Create ManyBodyCore
    core = ManyBodyCore(
        molecule=molecule,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf", 2: "mp2"},
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    print("✓ ManyBodyCore created successfully")

    # Check if dependency graph method is available
    if hasattr(core, 'iterate_molecules_by_level'):
        print("✓ iterate_molecules_by_level method available")
    else:
        print("✗ iterate_molecules_by_level method NOT available")
        return False

    # Create parallel configuration (testing mode without QCEngine)
    config = ParallelConfig(
        max_workers=2,
        execution_mode="threading",
        use_qcengine=False,  # Use placeholder calculations for testing
        timeout_seconds=60
    )

    print("✓ ParallelConfig created")

    specifications = {
        "hf": {
            "program": "psi4",
            "specification": {
                "driver": "energy",
                "model": {"method": "hf", "basis": "sto-3g"},
                "keywords": {},
                "protocols": {},
                "extras": {},
            },
        },
        "mp2": {
            "program": "psi4",
            "specification": {
                "driver": "energy",
                "model": {"method": "mp2", "basis": "sto-3g"},
                "keywords": {},
                "protocols": {},
                "extras": {},
            },
        },
    }

    # Create parallel executor
    try:
        executor = ParallelManyBodyExecutor(core, config, driver="energy", specifications=specifications)
        print("✓ ParallelManyBodyExecutor created successfully")
    except Exception as e:
        print(f"✗ Failed to create ParallelManyBodyExecutor: {e}")
        return False

    # Test parallel execution
    try:
        print("\nRunning parallel calculation...")
        results = executor.execute_full_calculation()
        print(f"✓ Parallel calculation completed: {len(results)} fragment results")

        # Display some results
        for label, result in list(results.items())[:3]:  # Show first 3 results
            print(f"  {label}: {result.return_result:.6f} Eh")

        # Show execution statistics
        stats = executor.get_execution_statistics()
        print(f"\nExecution Statistics:")
        print(f"  Total fragments: {stats['total_fragments']}")
        print(f"  Levels executed: {stats['levels_executed']}")
        print(f"  Parallel time: {stats['parallel_time']:.3f}s")
        print(f"  Estimated speedup: {stats['speedup_factor']:.2f}x")

        return True

    except Exception as e:
        print(f"✗ Parallel execution failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_sequential_comparison():
    """Test that parallel and sequential give identical results."""

    from qcmanybody import ManyBodyCore, BsseEnum
    from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
    from qcelemental.models import Molecule

    print("\nTesting Sequential vs Parallel Correctness...")
    print("="*50)

    # Simple dimer
    molecule = Molecule(**{
        'symbols': ['O', 'H', 'H', 'O', 'H', 'H'],
        'geometry': [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  1.43233673, -0.96104039],
            [ 0.0, -1.43233673, -0.96104039],
            [ 3.0,  0.0,  0.0],
            [ 3.0,  1.43233673, -0.96104039],
            [ 3.0, -1.43233673, -0.96104039],
        ],
        'fragments': [[0, 1, 2], [3, 4, 5]],
        'molecular_charge': 0.0,
        'molecular_multiplicity': 1
    })

    core = ManyBodyCore(
        molecule=molecule,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf", 2: "mp2"},
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    # Sequential execution
    sequential_config = ParallelConfig(
        max_workers=1,
        execution_mode="serial",
        use_qcengine=False
    )

    specifications = {
        "hf": {
            "program": "psi4",
            "specification": {
                "driver": "energy",
                "model": {"method": "hf", "basis": "sto-3g"},
                "keywords": {},
                "protocols": {},
                "extras": {},
            },
        },
        "mp2": {
            "program": "psi4",
            "specification": {
                "driver": "energy",
                "model": {"method": "mp2", "basis": "sto-3g"},
                "keywords": {},
                "protocols": {},
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
    print(f"✓ Sequential execution: {len(sequential_results)} results")

    # Parallel execution
    parallel_config = ParallelConfig(
        max_workers=2,
        execution_mode="threading",
        use_qcengine=False
    )

    parallel_executor = ParallelManyBodyExecutor(
        core,
        parallel_config,
        driver="energy",
        specifications=copy.deepcopy(specifications),
    )
    parallel_results = parallel_executor.execute_full_calculation()
    print(f"✓ Parallel execution: {len(parallel_results)} results")

    # Validate correctness
    try:
        parallel_executor.validate_parallel_correctness(
            parallel_results, sequential_results, tolerance=1e-12
        )
        print("✓ Parallel vs Sequential validation PASSED (tolerance=1e-12)")
        return True
    except Exception as e:
        print(f"✗ Parallel vs Sequential validation FAILED: {e}")
        return False


if __name__ == "__main__":
    print("QCManyBody Parallel Execution System Test")
    print("="*60)

    # Test 1: Basic functionality
    success1 = test_parallel_system()

    # Test 2: Sequential vs parallel correctness
    success2 = test_sequential_comparison()

    print("\n" + "="*60)
    if success1 and success2:
        print("✅ ALL TESTS PASSED - Parallel execution system is working!")
        print("\nThe parallel execution system is ready for use.")
        print("You can now run actual_test_water16.py for real QC calculations.")
    else:
        print("❌ SOME TESTS FAILED - Check the errors above.")
        print("\nThe parallel execution system needs fixes before use.")
    print("="*60)