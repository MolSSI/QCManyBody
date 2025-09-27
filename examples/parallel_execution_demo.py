#!/usr/bin/env python3
"""Parallel Execution Demo for QCManyBody.

This example demonstrates the use of the ParallelManyBodyExecutor for level-by-level
parallel execution of many-body calculations while respecting mathematical dependencies.
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
    """Demonstrate parallel execution capabilities."""
    print("🚀 QCManyBody Parallel Execution Demo")
    print("=" * 60)

    # Create a water dimer system
    print("\n📍 Creating water dimer system...")
    water_dimer = Molecule.from_data("""
    O  0.0000  0.0000  0.0000
    H  0.7570  0.5860  0.0000
    H -0.7570  0.5860  0.0000
    --
    O  3.0000  0.0000  0.0000
    H  3.7570  0.5860  0.0000
    H  2.2430  0.5860  0.0000
    """)

    # Create ManyBodyCore
    print("🧮 Setting up many-body calculation...")
    core = ManyBodyCore(
        molecule=water_dimer,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf", 2: "mp2"},  # Multi-level calculation
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    print(f"   - System: {len(water_dimer.fragments)} fragments")
    print(f"   - Levels: {core.levels}")
    print(f"   - BSSE treatment: {[bsse.value for bsse in core.bsse_type]}")

    # Configure parallel execution
    print("\n⚙️  Configuring parallel execution...")
    parallel_config = ParallelConfig(
        max_workers=4,
        execution_mode="threading",
        use_qcengine=False,  # Use placeholder for demo
        qc_program="psi4",
        basis_set="sto-3g"
    )

    print(f"   - Workers: {parallel_config.max_workers}")
    print(f"   - Mode: {parallel_config.execution_mode}")
    print(f"   - QCEngine: {'enabled' if parallel_config.use_qcengine else 'disabled (demo mode)'}")

    # Create parallel executor
    print("\n🔧 Initializing parallel executor...")
    executor = ParallelManyBodyExecutor(core, parallel_config)

    # Show dependency information
    print("\n📋 Dependency analysis:")
    print("   Level-by-level fragments:")
    level_count = {}
    for level, mc, label, mol in core.iterate_molecules_by_level():
        if level not in level_count:
            level_count[level] = 0
        level_count[level] += 1

    for level, count in level_count.items():
        print(f"     Level {level}: {count} fragments")

    # Execute parallel calculation
    print("\n⏱️  Running parallel calculation...")
    start_time = time.time()
    results = executor.execute_full_calculation()
    execution_time = time.time() - start_time

    # Show results
    print(f"\n✅ Parallel execution completed!")
    print(f"   - Execution time: {execution_time:.3f} seconds")
    print(f"   - Fragments calculated: {len(results)}")

    # Show execution statistics
    stats = executor.get_execution_statistics()
    print(f"   - Levels executed: {stats['levels_executed']}")
    print(f"   - Estimated speedup: {stats['speedup_factor']:.2f}x")

    # Show sample results
    print("\n📊 Sample results:")
    for i, (label, result) in enumerate(results.items()):
        if i < 3:  # Show first 3 results
            print(f"   {label}: Energy = {result.return_result:.6f}")
        elif i == 3:
            print(f"   ... and {len(results) - 3} more fragments")
            break

    # Demonstrate validation capabilities
    print("\n🔍 Validation capabilities:")
    print("   - Ultra-strict numerical precision (1e-12 tolerance)")
    print("   - Parallel vs sequential comparison")
    print("   - Mathematical dependency verification")
    print("   - Performance benchmarking")

    # Show configuration options
    print("\n🎛️  Available configuration options:")
    print("   - Execution modes: serial, threading, multiprocessing")
    print("   - QCEngine integration for real QC calculations")
    print("   - Flexible worker count and resource management")
    print("   - Error handling and recovery mechanisms")
    print("   - Performance monitoring and statistics")

    print("\n🎯 Demo completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()