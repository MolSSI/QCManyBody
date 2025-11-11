"""
Example demonstrating parallel execution in QCManyBody.

This example shows how to use the parallel execution module to
speed up many-body expansion calculations.
"""

import time
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput


def create_water_trimer():
    """Create a water trimer test system."""
    return Molecule(
        **{
            "symbols": ["O", "H", "H", "O", "H", "H", "O", "H", "H"],
            "geometry": [
                # Water 1
                [0.0, 0.0, 0.0],
                [0.0, 1.5, 0.0],
                [1.5, 0.0, 0.0],
                # Water 2
                [5.0, 0.0, 0.0],
                [5.0, 1.5, 0.0],
                [6.5, 0.0, 0.0],
                # Water 3
                [2.5, 4.0, 0.0],
                [2.5, 5.5, 0.0],
                [4.0, 4.0, 0.0],
            ],
            "fragments": [[0, 1, 2], [3, 4, 5], [6, 7, 8]],
        }
    )


def run_sequential(mol):
    """Run calculation sequentially (no parallelism)."""
    print("=" * 60)
    print("Running SEQUENTIAL calculation...")
    print("=" * 60)

    mbin = ManyBodyInput(
        molecule=mol,
        specification={
            "driver": "energy",
            "keywords": {"bsse_type": ["cp"], "max_nbody": 2},
            "specification": {
                "hf/sto-3g": {"program": "psi4", "model": {"method": "hf", "basis": "sto-3g"}}
            },
        },
    )

    start = time.time()
    result = ManyBodyComputer.from_manybodyinput(mbin, parallel=False)
    elapsed = time.time() - start

    print(f"\n✓ Sequential calculation complete!")
    print(f"  Time: {elapsed:.2f}s")
    print(f"  Energy: {result.return_result:.8f} Eh")

    return elapsed, result


def run_parallel_basic(mol):
    """Run calculation with basic parallel settings."""
    print("\n" + "=" * 60)
    print("Running PARALLEL calculation (auto-detect cores)...")
    print("=" * 60)

    mbin = ManyBodyInput(
        molecule=mol,
        specification={
            "driver": "energy",
            "keywords": {"bsse_type": ["cp"], "max_nbody": 2},
            "specification": {
                "hf/sto-3g": {"program": "psi4", "model": {"method": "hf", "basis": "sto-3g"}}
            },
        },
    )

    start = time.time()
    result = ManyBodyComputer.from_manybodyinput(mbin, parallel=True)
    elapsed = time.time() - start

    print(f"\n✓ Parallel calculation complete!")
    print(f"  Time: {elapsed:.2f}s")
    print(f"  Energy: {result.return_result:.8f} Eh")

    return elapsed, result


def run_parallel_custom(mol, n_workers=4):
    """Run calculation with custom parallel configuration."""
    print("\n" + "=" * 60)
    print(f"Running PARALLEL calculation ({n_workers} workers)...")
    print("=" * 60)

    from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

    # Custom configuration
    config = ExecutorConfig(
        n_workers=n_workers,
        timeout_per_task=300,  # 5 minutes per task
        max_retries=2,
        log_level="INFO",
    )

    executor = MultiprocessingExecutor(config)

    mbin = ManyBodyInput(
        molecule=mol,
        specification={
            "driver": "energy",
            "keywords": {"bsse_type": ["cp"], "max_nbody": 2},
            "specification": {
                "hf/sto-3g": {"program": "psi4", "model": {"method": "hf", "basis": "sto-3g"}}
            },
        },
    )

    start = time.time()

    with executor:
        print(f"Executor info: {executor.get_info()}")
        result = ManyBodyComputer.from_manybodyinput(mbin, executor=executor)

    elapsed = time.time() - start

    print(f"\n✓ Parallel calculation complete!")
    print(f"  Time: {elapsed:.2f}s")
    print(f"  Energy: {result.return_result:.8f} Eh")

    return elapsed, result


def main():
    """Run all examples and compare performance."""
    print("\n" + "=" * 60)
    print("QCManyBody Parallel Execution Example")
    print("=" * 60)

    # Create test molecule
    mol = create_water_trimer()
    print(f"\nTest system: Water trimer")
    print(f"  Fragments: {len(mol.fragments)}")
    print(f"  Atoms: {len(mol.symbols)}")

    try:
        # Run sequential
        time_seq, result_seq = run_sequential(mol)

        # Run parallel (basic)
        time_par_basic, result_par_basic = run_parallel_basic(mol)

        # Run parallel (custom)
        time_par_custom, result_par_custom = run_parallel_custom(mol, n_workers=2)

        # Compare results
        print("\n" + "=" * 60)
        print("PERFORMANCE COMPARISON")
        print("=" * 60)
        print(f"Sequential:          {time_seq:.2f}s")
        print(f"Parallel (auto):     {time_par_basic:.2f}s  (speedup: {time_seq/time_par_basic:.2f}x)")
        print(f"Parallel (2 workers): {time_par_custom:.2f}s  (speedup: {time_seq/time_par_custom:.2f}x)")

        # Verify results match
        print("\n" + "=" * 60)
        print("RESULT VERIFICATION")
        print("=" * 60)
        print(f"Sequential energy:        {result_seq.return_result:.8f} Eh")
        print(f"Parallel (auto) energy:   {result_par_basic.return_result:.8f} Eh")
        print(f"Parallel (custom) energy: {result_par_custom.return_result:.8f} Eh")

        diff1 = abs(result_seq.return_result - result_par_basic.return_result)
        diff2 = abs(result_seq.return_result - result_par_custom.return_result)

        print(f"\nDifference (seq vs par_basic):  {diff1:.2e}")
        print(f"Difference (seq vs par_custom): {diff2:.2e}")

        if diff1 < 1e-10 and diff2 < 1e-10:
            print("\n✓ Results match within tolerance!")
        else:
            print("\n⚠ Results differ more than expected")

    except ImportError as e:
        print(f"\n⚠ Error: {e}")
        print("This example requires qcengine and psi4 to be installed.")
        print("Install with: conda install -c conda-forge qcengine psi4")
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
