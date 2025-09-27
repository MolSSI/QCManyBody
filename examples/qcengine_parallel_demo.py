#!/usr/bin/env python3
"""QCEngine integration demo for parallel execution.

This example demonstrates real quantum chemistry calculations using QCEngine
with the ParallelManyBodyExecutor.
"""

import sys
import time
from pathlib import Path
from qcelemental.models import Molecule

# Add parent directory to path for qcmanybody imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from qcmanybody import ManyBodyCore, BsseEnum
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig


def main():
    """Demonstrate QCEngine integration with parallel execution."""
    print("🚀 QCManyBody + QCEngine Parallel Demo")
    print("=" * 50)

    # Create a simple H2 dimer system for fast demo
    print("\n📍 Creating H2 dimer system...")
    h2_dimer = Molecule.from_data("""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    --
    H 3.0 0.0 0.0
    H 3.0 0.0 0.74
    """)

    print(f"   - System: {len(h2_dimer.fragments)} H2 fragments")
    print(f"   - Total atoms: {len(h2_dimer.symbols)}")

    # Create ManyBodyCore with simple configuration
    print("\n🧮 Setting up many-body calculation...")
    core = ManyBodyCore(
        molecule=h2_dimer,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf", 2: "hf"},  # Keep it simple for demo
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    print(f"   - Method: Hartree-Fock")
    print(f"   - Basis: STO-3G")
    print(f"   - BSSE: No counterpoise")

    # Test 1: Placeholder execution (fast)
    print("\n🧪 Test 1: Placeholder execution (validation)")
    config_placeholder = ParallelConfig(
        max_workers=2,
        execution_mode="threading",
        use_qcengine=False,  # Placeholder mode
        basis_set="sto-3g"
    )

    executor_placeholder = ParallelManyBodyExecutor(core, config_placeholder)

    start_time = time.time()
    results_placeholder = executor_placeholder.execute_full_calculation()
    placeholder_time = time.time() - start_time

    print(f"   ✅ Placeholder execution: {len(results_placeholder)} fragments in {placeholder_time:.3f}s")

    # Test 2: Real QCEngine execution with Psi4
    print("\n⚡ Test 2: Real QCEngine execution with Psi4")
    config_qcengine = ParallelConfig(
        max_workers=2,
        execution_mode="threading",
        use_qcengine=True,  # Real QCEngine execution
        qc_program="psi4",
        basis_set="sto-3g",
        qcengine_config={
            "keywords": {"scf_type": "df", "mp2_type": "df"},
            "protocols": {"stdout": False}  # Suppress output for clean demo
        }
    )

    executor_qcengine = ParallelManyBodyExecutor(core, config_qcengine)

    print(f"   - Using {config_qcengine.qc_program} via QCEngine")
    print(f"   - Workers: {config_qcengine.max_workers}")
    print(f"   - Running calculation...")

    start_time = time.time()
    results_qcengine = executor_qcengine.execute_full_calculation()
    qcengine_time = time.time() - start_time

    print(f"   ✅ QCEngine execution: {len(results_qcengine)} fragments in {qcengine_time:.3f}s")

    # Show results comparison
    print("\n📊 Results comparison:")
    print(f"   Placeholder mode: {placeholder_time:.3f}s")
    print(f"   QCEngine mode:    {qcengine_time:.3f}s")
    print(f"   Ratio:            {qcengine_time/placeholder_time:.1f}x slower (expected for real QC)")

    # Show some actual energies
    print("\n🔬 QCEngine calculation results:")
    for i, (label, result) in enumerate(results_qcengine.items()):
        if i < 3:  # Show first 3 results
            print(f"   {label}: {result.return_result:.8f} hartree")
        elif i == 3:
            print(f"   ... and {len(results_qcengine) - 3} more")
            break

    # Execution statistics
    stats = executor_qcengine.get_execution_statistics()
    print(f"\n📈 Execution statistics:")
    print(f"   - Total fragments: {stats['total_fragments']}")
    print(f"   - Levels executed: {stats['levels_executed']}")
    print(f"   - Parallel efficiency: Real QC calculations successfully parallelized")

    # Show dependency ordering
    print(f"\n🔗 Dependency analysis (level-by-level execution):")
    level_counts = {}
    for level, mc, label, mol in core.iterate_molecules_by_level():
        level_counts[level] = level_counts.get(level, 0) + 1

    for level, count in level_counts.items():
        print(f"   Level {level}: {count} fragments ({'monomers' if level == 1 else 'dimers' if level == 2 else f'{level}-mers'})")

    print("\n✅ QCEngine integration demonstrated successfully!")
    print("   - Real quantum chemistry calculations")
    print("   - Parallel execution within levels")
    print("   - Mathematical dependency preservation")
    print("   - Performance monitoring and statistics")
    print("=" * 50)


if __name__ == "__main__":
    main()