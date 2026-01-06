#!/usr/bin/env python3
"""Parallel execution performance benchmarking tool.

This script provides comprehensive performance analysis of the ParallelManyBodyExecutor,
including scalability testing, speedup analysis, and validation against sequential execution.
"""

import argparse
import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from qcelemental.models import Molecule

# Add parent directory to path for qcmanybody imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from qcmanybody.core import ManyBodyCore
from qcmanybody.parallel import ParallelManyBodyExecutor, ParallelConfig
from qcmanybody.models.v1 import BsseEnum


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def create_test_systems() -> Dict[str, Molecule]:
    """Create a set of test systems for benchmarking.

    Returns
    -------
    Dict[str, Molecule]
        Dictionary mapping system names to molecules
    """
    systems = {}

    # Small system: Water dimer
    systems["water_dimer"] = Molecule.from_data("""
    O  0.0000  0.0000  0.0000
    H  0.7570  0.5860  0.0000
    H -0.7570  0.5860  0.0000
    --
    O  3.0000  0.0000  0.0000
    H  3.7570  0.5860  0.0000
    H  2.2430  0.5860  0.0000
    """)

    # Medium system: Water trimer
    systems["water_trimer"] = Molecule.from_data("""
    O  0.0000  0.0000  0.0000
    H  0.7570  0.5860  0.0000
    H -0.7570  0.5860  0.0000
    --
    O  3.0000  0.0000  0.0000
    H  3.7570  0.5860  0.0000
    H  2.2430  0.5860  0.0000
    --
    O  0.0000  3.0000  0.0000
    H  0.7570  3.5860  0.0000
    H -0.7570  3.5860  0.0000
    """)

    # Large system: Water tetramer
    systems["water_tetramer"] = Molecule.from_data("""
    O  0.0000  0.0000  0.0000
    H  0.7570  0.5860  0.0000
    H -0.7570  0.5860  0.0000
    --
    O  3.0000  0.0000  0.0000
    H  3.7570  0.5860  0.0000
    H  2.2430  0.5860  0.0000
    --
    O  0.0000  3.0000  0.0000
    H  0.7570  3.5860  0.0000
    H -0.7570  3.5860  0.0000
    --
    O  3.0000  3.0000  0.0000
    H  3.7570  3.5860  0.0000
    H  2.2430  3.5860  0.0000
    """)

    return systems


def create_manybody_core(molecule: Molecule, max_nbody: int = 2) -> ManyBodyCore:
    """Create a ManyBodyCore for benchmarking.

    Parameters
    ----------
    molecule : Molecule
        The molecular system
    max_nbody : int
        Maximum N-body level to compute

    Returns
    -------
    ManyBodyCore
        Configured ManyBodyCore instance
    """
    nfragments = len(molecule.fragments)
    max_nbody = min(max_nbody, nfragments)

    levels = {i: "hf" for i in range(1, max_nbody + 1)}

    return ManyBodyCore(
        molecule=molecule,
        bsse_type=[BsseEnum.nocp],
        levels=levels,
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={}
    )


def benchmark_sequential_execution(core: ManyBodyCore) -> Tuple[Dict, float]:
    """Benchmark sequential execution using original iterate_molecules.

    Parameters
    ----------
    core : ManyBodyCore
        ManyBodyCore instance to benchmark

    Returns
    -------
    Tuple[Dict, float]
        Sequential results and execution time
    """
    start_time = time.time()

    # Simulate sequential execution by collecting all fragments
    sequential_results = {}
    fragment_count = 0

    for mc, label, mol in core.iterate_molecules():
        fragment_count += 1
        # Simulate calculation time
        time.sleep(0.01)

        # Create placeholder result
        natoms = len(mol.symbols)
        energy_estimate = -natoms * 1.0

        sequential_results[label] = {
            "energy": energy_estimate,
            "mc": mc,
            "natoms": natoms
        }

    execution_time = time.time() - start_time

    logging.info(f"Sequential execution: {fragment_count} fragments in {execution_time:.3f}s")

    return sequential_results, execution_time


def benchmark_parallel_execution(core: ManyBodyCore, config: ParallelConfig) -> Tuple[Dict, Dict]:
    """Benchmark parallel execution.

    Parameters
    ----------
    core : ManyBodyCore
        ManyBodyCore instance to benchmark
    config : ParallelConfig
        Parallel execution configuration

    Returns
    -------
    Tuple[Dict, Dict]
        Parallel results and execution statistics
    """
    executor = ParallelManyBodyExecutor(
        core,
        config,
        driver=config.default_driver,
        specifications={},
    )

    start_time = time.time()
    parallel_results = executor.execute_full_calculation()
    execution_time = time.time() - start_time

    stats = executor.get_execution_statistics()
    stats["actual_execution_time"] = execution_time

    logging.info(f"Parallel execution ({config.execution_mode}, {config.max_workers} workers): "
                f"{stats['total_fragments']} fragments in {execution_time:.3f}s")

    return parallel_results, stats


def validate_correctness(sequential_results: Dict, parallel_results: Dict, tolerance: float = 1e-12) -> bool:
    """Validate parallel results against sequential results.

    Parameters
    ----------
    sequential_results : Dict
        Results from sequential execution
    parallel_results : Dict
        Results from parallel execution
    tolerance : float
        Numerical tolerance for comparison

    Returns
    -------
    bool
        True if results match within tolerance
    """
    # Convert parallel AtomicResult objects to comparable format
    parallel_simplified = {}
    for label, result in parallel_results.items():
        model_method = getattr(result.model, "method", None)
        if model_method is None and isinstance(result.model, dict):
            model_method = result.model.get("method")
        parallel_simplified[label] = {
            "energy": result.return_result,
            "mc": model_method,
            "natoms": len(result.molecule.symbols)
        }

    # Check label sets match
    seq_labels = set(sequential_results.keys())
    par_labels = set(parallel_simplified.keys())

    if seq_labels != par_labels:
        logging.error(f"Label mismatch: sequential={len(seq_labels)}, parallel={len(par_labels)}")
        return False

    # Check numerical values
    max_difference = 0.0
    for label in seq_labels:
        seq_energy = sequential_results[label]["energy"]
        par_energy = parallel_simplified[label]["energy"]
        difference = abs(seq_energy - par_energy)
        max_difference = max(max_difference, difference)

        if difference > tolerance:
            logging.error(f"Energy difference for {label}: {difference} > {tolerance}")
            return False

    logging.info(f"Validation passed: max difference = {max_difference}")
    return True


def analyze_scalability(system_name: str, molecule: Molecule, worker_counts: List[int]) -> Dict:
    """Analyze scalability across different worker counts.

    Parameters
    ----------
    system_name : str
        Name of the test system
    molecule : Molecule
        Molecular system to test
    worker_counts : List[int]
        List of worker counts to test

    Returns
    -------
    Dict
        Scalability analysis results
    """
    logging.info(f"Analyzing scalability for {system_name}")

    core = create_manybody_core(molecule)

    # Benchmark sequential execution first
    sequential_results, sequential_time = benchmark_sequential_execution(core)

    scalability_results = {
        "system": system_name,
        "nfragments": len(molecule.fragments),
        "sequential_time": sequential_time,
        "sequential_fragment_count": len(sequential_results),
        "worker_results": {}
    }

    for worker_count in worker_counts:
        logging.info(f"Testing with {worker_count} workers")

        config = ParallelConfig(
            max_workers=worker_count,
            execution_mode="threading",  # Use threading for consistent comparison
            use_qcengine=False  # Use placeholder for consistent timing
        )

        parallel_results, stats = benchmark_parallel_execution(core, config)

        # Validate correctness
        is_correct = validate_correctness(sequential_results, parallel_results)

        # Calculate speedup
        if stats.get("speedup_factor") is not None:
            speedup = stats["speedup_factor"]
        elif stats["actual_execution_time"] > 0:
            speedup = sequential_time / stats["actual_execution_time"]
        else:
            speedup = 0

        # Calculate efficiency
        efficiency = speedup / worker_count if worker_count > 0 else 0

        serial_time = stats.get("serial_time")
        if serial_time is None and stats.get("speedup_factor") is not None:
            serial_time = stats["speedup_factor"] * stats["actual_execution_time"]
        if serial_time is None:
            serial_time = sequential_time

        scalability_results["worker_results"][worker_count] = {
            "execution_time": stats["actual_execution_time"],
            "fragment_count": stats["total_fragments"],
            "levels_executed": stats["levels_executed"],
            "speedup": speedup,
            "efficiency": efficiency,
            "is_correct": is_correct,
            "serial_time": serial_time,
            "parallel_overhead": serial_time - stats["actual_execution_time"]
        }

        logging.info(f"  Workers: {worker_count}, Time: {stats['actual_execution_time']:.3f}s, "
                    f"Speedup: {speedup:.2f}x, Efficiency: {efficiency:.2f}")

    return scalability_results


def run_comprehensive_benchmark(output_file: str = "parallel_benchmark_results.json"):
    """Run comprehensive performance benchmark.

    Parameters
    ----------
    output_file : str
        Path to save benchmark results
    """
    logging.info("Starting comprehensive parallel execution benchmark")

    systems = create_test_systems()
    worker_counts = [1, 2, 4, 8]

    all_results = {
        "benchmark_timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "systems_tested": list(systems.keys()),
        "worker_counts_tested": worker_counts,
        "results": {}
    }

    for system_name, molecule in systems.items():
        logging.info(f"\n{'='*60}")
        logging.info(f"Benchmarking system: {system_name}")
        logging.info(f"Fragments: {len(molecule.fragments)}")

        try:
            scalability_results = analyze_scalability(system_name, molecule, worker_counts)
            all_results["results"][system_name] = scalability_results

            # Print summary for this system
            best_speedup = max(
                result["speedup"] for result in scalability_results["worker_results"].values()
            )
            best_workers = max(
                scalability_results["worker_results"].items(),
                key=lambda x: x[1]["speedup"]
            )[0]

            logging.info(f"Best speedup for {system_name}: {best_speedup:.2f}x with {best_workers} workers")

        except Exception as e:
            logging.error(f"Benchmark failed for {system_name}: {e}")
            all_results["results"][system_name] = {"error": str(e)}

    # Save results
    output_path = Path(output_file)
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)

    logging.info(f"\nBenchmark results saved to: {output_path}")

    # Print summary
    print_benchmark_summary(all_results)


def print_benchmark_summary(results: Dict):
    """Print a summary of benchmark results.

    Parameters
    ----------
    results : Dict
        Benchmark results to summarize
    """
    print(f"\n{'='*80}")
    print("PARALLEL EXECUTION BENCHMARK SUMMARY")
    print(f"{'='*80}")

    for system_name, system_results in results["results"].items():
        if "error" in system_results:
            print(f"\n{system_name}: FAILED - {system_results['error']}")
            continue

        print(f"\n{system_name.upper()}:")
        print(f"  Fragments: {system_results['nfragments']}")
        print(f"  Sequential time: {system_results['sequential_time']:.3f}s")

        print(f"  {'Workers':<8} {'Time (s)':<10} {'Speedup':<10} {'Efficiency':<12} {'Correct':<8}")
        print(f"  {'-'*48}")

        for workers, worker_result in system_results["worker_results"].items():
            print(f"  {workers:<8} {worker_result['execution_time']:<10.3f} "
                  f"{worker_result['speedup']:<10.2f} {worker_result['efficiency']:<12.2f} "
                  f"{'✓' if worker_result['is_correct'] else '✗':<8}")

    print(f"\n{'='*80}")


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description="Benchmark parallel execution performance")
    parser.add_argument("--output", default="parallel_benchmark_results.json",
                       help="Output file for benchmark results")
    parser.add_argument("--systems", nargs="+",
                       choices=["water_dimer", "water_trimer", "water_tetramer", "all"],
                       default=["all"],
                       help="Systems to benchmark")
    parser.add_argument("--workers", nargs="+", type=int, default=[1, 2, 4, 8],
                       help="Worker counts to test")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose logging")

    args = parser.parse_args()

    setup_logging(args.verbose)

    if "all" in args.systems:
        run_comprehensive_benchmark(args.output)
    else:
        # Custom benchmark with selected systems
        systems = create_test_systems()
        selected_systems = {name: systems[name] for name in args.systems if name in systems}

        all_results = {
            "benchmark_timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "systems_tested": list(selected_systems.keys()),
            "worker_counts_tested": args.workers,
            "results": {}
        }

        for system_name, molecule in selected_systems.items():
            scalability_results = analyze_scalability(system_name, molecule, args.workers)
            all_results["results"][system_name] = scalability_results

        # Save and print results
        with open(args.output, 'w') as f:
            json.dump(all_results, f, indent=2)

        print_benchmark_summary(all_results)


if __name__ == "__main__":
    main()