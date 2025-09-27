#!/usr/bin/env python3
"""
Performance Benchmarking for QCManyBody Dependency Graph Implementation

This script comprehensively benchmarks the dependency graph performance
for various system sizes and configurations, providing detailed analysis
for P1-002 performance optimization requirements.
"""

import time
import sys
import tracemalloc
from typing import Dict, List, Tuple, Any
import psutil
import gc
import qcelemental as qcel
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum


def create_large_molecule(n_fragments: int) -> qcel.models.Molecule:
    """Create a test molecule with specified number of fragments.

    Parameters
    ----------
    n_fragments : int
        Number of fragments to create

    Returns
    -------
    qcel.models.Molecule
        Test molecule with n_fragments
    """
    fragments = []
    for i in range(n_fragments):
        # Simple helium atoms separated in space
        x = i * 5.0  # 5 Bohr separation
        fragments.append(f"He {x} 0 0")

    molecule_string = "\n--\n".join(fragments)
    return qcel.models.Molecule.from_data(molecule_string)


def benchmark_dependency_graph_construction(n_fragments: int,
                                          max_nbody: int = 2,
                                          n_runs: int = 5) -> Dict[str, Any]:
    """Benchmark dependency graph construction time and memory usage.

    Parameters
    ----------
    n_fragments : int
        Number of fragments in test system
    max_nbody : int
        Maximum N-body level
    n_runs : int
        Number of benchmark runs for averaging

    Returns
    -------
    Dict[str, Any]
        Benchmark results
    """
    molecule = create_large_molecule(n_fragments)

    # Setup configuration
    levels = {i: "hf" for i in range(1, max_nbody + 1)}

    construction_times = []
    memory_peaks = []

    for run in range(n_runs):
        # Start memory tracking
        tracemalloc.start()
        gc.collect()  # Clean slate

        # Benchmark construction
        start_time = time.perf_counter()

        mbc = ManyBodyCore(
            molecule=molecule,
            bsse_type=[BsseEnum.cp],
            levels=levels,
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Access dependency graph (triggers construction)
        dep_graph = mbc.dependency_graph

        end_time = time.perf_counter()

        # Measure memory
        current, peak = tracemalloc.get_traced_memory()
        memory_peaks.append(peak / 1024 / 1024)  # Convert to MB
        tracemalloc.stop()

        construction_times.append(end_time - start_time)

        # Cleanup
        del mbc, dep_graph
        gc.collect()

    return {
        "n_fragments": n_fragments,
        "max_nbody": max_nbody,
        "construction_time_avg": sum(construction_times) / len(construction_times),
        "construction_time_std": (sum((t - sum(construction_times)/len(construction_times))**2 for t in construction_times) / len(construction_times))**0.5,
        "memory_peak_avg_mb": sum(memory_peaks) / len(memory_peaks),
        "memory_peak_std_mb": (sum((m - sum(memory_peaks)/len(memory_peaks))**2 for m in memory_peaks) / len(memory_peaks))**0.5,
        "total_fragments_analyzed": len(list(mbc.iterate_molecules())) if 'mbc' in locals() else 0
    }


def benchmark_iteration_performance(n_fragments: int,
                                   max_nbody: int = 2,
                                   n_runs: int = 10) -> Dict[str, Any]:
    """Benchmark iteration performance comparison.

    Parameters
    ----------
    n_fragments : int
        Number of fragments in test system
    max_nbody : int
        Maximum N-body level
    n_runs : int
        Number of benchmark runs

    Returns
    -------
    Dict[str, Any]
        Performance comparison results
    """
    molecule = create_large_molecule(n_fragments)
    levels = {i: "hf" for i in range(1, max_nbody + 1)}

    mbc = ManyBodyCore(
        molecule=molecule,
        bsse_type=[BsseEnum.cp],
        levels=levels,
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    # Benchmark original iterate_molecules()
    original_times = []
    for run in range(n_runs):
        start_time = time.perf_counter()
        fragments_original = list(mbc.iterate_molecules())
        end_time = time.perf_counter()
        original_times.append(end_time - start_time)

    # Benchmark new iterate_molecules_by_level()
    level_times = []
    for run in range(n_runs):
        start_time = time.perf_counter()
        fragments_level = list(mbc.iterate_molecules_by_level())
        end_time = time.perf_counter()
        level_times.append(end_time - start_time)

    return {
        "n_fragments": n_fragments,
        "max_nbody": max_nbody,
        "original_method_time_avg": sum(original_times) / len(original_times),
        "level_method_time_avg": sum(level_times) / len(level_times),
        "overhead_percentage": ((sum(level_times) / len(level_times)) - (sum(original_times) / len(original_times))) / (sum(original_times) / len(original_times)) * 100,
        "fragments_count": len(fragments_original),
        "level_fragments_count": len(fragments_level)
    }


def analyze_scalability(max_fragments: int = 20,
                       step_size: int = 2) -> List[Dict[str, Any]]:
    """Analyze scalability across different system sizes.

    Parameters
    ----------
    max_fragments : int
        Maximum number of fragments to test
    step_size : int
        Step size for fragment count

    Returns
    -------
    List[Dict[str, Any]]
        Scalability analysis results
    """
    results = []

    print(f"🔍 Analyzing scalability from 2 to {max_fragments} fragments...")

    for n_frags in range(2, max_fragments + 1, step_size):
        print(f"  Testing {n_frags} fragments...", end="", flush=True)

        try:
            # Test with 2-body calculations
            construction_result = benchmark_dependency_graph_construction(n_frags, max_nbody=2, n_runs=3)
            iteration_result = benchmark_iteration_performance(n_frags, max_nbody=2, n_runs=5)

            combined_result = {**construction_result, **iteration_result}
            results.append(combined_result)

            print(f" ✓ (Construction: {construction_result['construction_time_avg']:.4f}s, Overhead: {iteration_result['overhead_percentage']:.1f}%)")

        except Exception as e:
            print(f" ❌ Error: {e}")
            continue

    return results


def benchmark_memory_usage(n_fragments_list: List[int]) -> List[Dict[str, Any]]:
    """Detailed memory usage analysis.

    Parameters
    ----------
    n_fragments_list : List[int]
        List of fragment counts to test

    Returns
    -------
    List[Dict[str, Any]]
        Memory usage analysis results
    """
    results = []

    print("\n💾 Analyzing memory usage patterns...")

    for n_frags in n_fragments_list:
        print(f"  Memory analysis for {n_frags} fragments...", end="", flush=True)

        # Get system memory before
        process = psutil.Process()
        memory_before = process.memory_info().rss / 1024 / 1024  # MB

        # Benchmark with detailed memory tracking
        tracemalloc.start()

        molecule = create_large_molecule(n_frags)
        levels = {1: "hf", 2: "mp2"}

        mbc = ManyBodyCore(
            molecule=molecule,
            bsse_type=[BsseEnum.cp],
            levels=levels,
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Create dependency graph
        dep_graph = mbc.dependency_graph

        # Iterate through all fragments to get peak usage
        fragments = list(mbc.iterate_molecules_by_level())

        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        memory_after = process.memory_info().rss / 1024 / 1024  # MB

        result = {
            "n_fragments": n_frags,
            "memory_peak_tracemalloc_mb": peak / 1024 / 1024,
            "memory_before_mb": memory_before,
            "memory_after_mb": memory_after,
            "memory_delta_mb": memory_after - memory_before,
            "fragments_processed": len(fragments),
            "dependency_levels": len(dep_graph.get_dependency_levels())
        }

        results.append(result)
        print(f" ✓ Peak: {result['memory_peak_tracemalloc_mb']:.1f}MB, Delta: {result['memory_delta_mb']:.1f}MB")

        # Cleanup
        del mbc, dep_graph, fragments, molecule
        gc.collect()

    return results


def generate_performance_report(scalability_results: List[Dict[str, Any]],
                               memory_results: List[Dict[str, Any]]) -> str:
    """Generate comprehensive performance report.

    Parameters
    ----------
    scalability_results : List[Dict[str, Any]]
        Results from scalability analysis
    memory_results : List[Dict[str, Any]]
        Results from memory analysis

    Returns
    -------
    str
        Formatted performance report
    """
    report = []
    report.append("="*80)
    report.append("🚀 QCManyBody Dependency Graph Performance Analysis Report")
    report.append("="*80)

    # Executive Summary
    if scalability_results:
        max_overhead = max(r['overhead_percentage'] for r in scalability_results)
        avg_overhead = sum(r['overhead_percentage'] for r in scalability_results) / len(scalability_results)

        report.append(f"\n📊 Executive Summary:")
        report.append(f"   • Systems tested: {len(scalability_results)} different fragment counts")
        report.append(f"   • Maximum performance overhead: {max_overhead:.1f}%")
        report.append(f"   • Average performance overhead: {avg_overhead:.1f}%")
        report.append(f"   • Largest system tested: {max(r['n_fragments'] for r in scalability_results)} fragments")

        # Performance acceptance
        if max_overhead < 5.0:
            report.append(f"   ✅ PERFORMANCE TARGET MET: Overhead < 5% requirement satisfied")
        else:
            report.append(f"   ⚠️  PERFORMANCE ATTENTION: Overhead {max_overhead:.1f}% exceeds 5% target")

    # Scalability Analysis
    report.append(f"\n🔍 Scalability Analysis:")
    report.append(f"{'Fragments':<10} {'Construction(s)':<15} {'Overhead(%)':<12} {'Total Fragments':<15}")
    report.append("-" * 60)

    for result in scalability_results:
        report.append(f"{result['n_fragments']:<10} "
                     f"{result['construction_time_avg']:<15.4f} "
                     f"{result['overhead_percentage']:<12.1f} "
                     f"{result['fragments_count']:<15}")

    # Memory Analysis
    if memory_results:
        report.append(f"\n💾 Memory Usage Analysis:")
        report.append(f"{'Fragments':<10} {'Peak (MB)':<12} {'Delta (MB)':<12} {'Levels':<8}")
        report.append("-" * 50)

        for result in memory_results:
            report.append(f"{result['n_fragments']:<10} "
                         f"{result['memory_peak_tracemalloc_mb']:<12.1f} "
                         f"{result['memory_delta_mb']:<12.1f} "
                         f"{result['dependency_levels']:<8}")

    # Recommendations
    report.append(f"\n🎯 Performance Optimization Recommendations:")

    if scalability_results:
        if max_overhead > 5.0:
            report.append(f"   • Optimize dependency graph construction for {max_overhead:.1f}% overhead")
            report.append(f"   • Consider lazy evaluation for large fragment systems")
            report.append(f"   • Profile memory allocation patterns in critical paths")
        else:
            report.append(f"   ✅ Current implementation meets performance targets")
            report.append(f"   • Consider further optimization for very large systems (50+ fragments)")

    if memory_results:
        max_memory = max(r['memory_peak_tracemalloc_mb'] for r in memory_results)
        report.append(f"   • Maximum memory usage: {max_memory:.1f}MB")
        if max_memory > 100:
            report.append(f"   • Consider memory optimization for systems requiring >100MB")

    # Phase 1 Task P1-002 Status
    report.append(f"\n🏁 Phase 1 Task P1-002 Status:")
    if scalability_results and max_overhead < 5.0:
        report.append(f"   ✅ Performance optimization requirement: SATISFIED")
    else:
        report.append(f"   ⚠️  Performance optimization requirement: NEEDS ATTENTION")

    report.append(f"   ✅ Scalability analysis: COMPLETE")
    report.append(f"   ✅ Memory profiling: COMPLETE")
    report.append(f"   ✅ Performance benchmarks: COMPLETE")

    report.append("="*80)

    return "\n".join(report)


def main():
    """Main benchmarking execution."""
    print("🚀 QCManyBody Dependency Graph Performance Benchmarking")
    print("=" * 60)

    # Run scalability analysis
    scalability_results = analyze_scalability(max_fragments=16, step_size=2)

    # Run memory analysis
    memory_test_sizes = [2, 4, 8, 12, 16, 20]
    memory_results = benchmark_memory_usage(memory_test_sizes)

    # Generate and display report
    report = generate_performance_report(scalability_results, memory_results)
    print(report)

    # Save detailed results
    import json

    results_data = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "scalability_results": scalability_results,
        "memory_results": memory_results,
        "system_info": {
            "python_version": sys.version,
            "total_memory_gb": psutil.virtual_memory().total / 1024**3,
            "cpu_count": psutil.cpu_count()
        }
    }

    with open("dependency_graph_performance_report.json", "w") as f:
        json.dump(results_data, f, indent=2)

    print(f"\n📄 Detailed results saved to: dependency_graph_performance_report.json")

    # Return success/failure for P1-002 completion
    if scalability_results:
        max_overhead = max(r['overhead_percentage'] for r in scalability_results)
        if max_overhead < 5.0:
            print(f"\n🎉 P1-002 PERFORMANCE REQUIREMENTS: ✅ SATISFIED")
            return 0
        else:
            print(f"\n⚠️  P1-002 PERFORMANCE REQUIREMENTS: ⚠️  NEEDS OPTIMIZATION")
            return 1
    else:
        print(f"\n❌ P1-002 PERFORMANCE REQUIREMENTS: ❌ BENCHMARKING FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())