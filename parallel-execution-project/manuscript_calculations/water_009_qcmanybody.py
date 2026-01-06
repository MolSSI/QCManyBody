#!/usr/bin/env python3
from __future__ import annotations

import time
from typing import Dict

from qcelemental.models import Molecule
from qcmanybody.models import (
    AtomicSpecification,
    ManyBodyInput,
    ManyBodyKeywords,
    ManyBodySpecification,
)
from qcmanybody.parallel import ParallelConfig, ParallelManyBodyExecutor

SYMBOLS = (
    ['O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H',
     'O',
     'H',
     'H']
)
GEOMETRY = (
    [[3.576212, 3.557473, 4.299181],
     [3.893906, 2.682676, 4.664189],
     [3.680219, 4.182103, 5.069044],
     [-0.01481, 3.032754, 2.640785],
     [-0.205779, 3.822767, 2.108249],
     [-0.779399, 2.484361, 2.456379],
     [2.234298, 1.342624, 1.324981],
     [1.896033, 1.926584, 0.61375],
     [2.841685, 1.98883, 1.701723],
     [-1.336477, 0.574082, 2.431147],
     [-0.412391, 0.325897, 2.270665],
     [-1.878907, -0.134295, 2.031022],
     [-4.127138, -3.832039, 2.076922],
     [-4.077673, -2.920405, 2.20451],
     [-4.420585, -4.178142, 2.934112],
     [-2.459635, -1.130859, -0.758204],
     [-2.722983, -1.994357, -0.991546],
     [-1.58271, -1.261066, -0.539542],
     [-4.161542, -0.786949, 1.708313],
     [-4.968122, -0.21796, 1.664925],
     [-3.648582, -0.404935, 1.039042],
     [-1.764704, -4.156586, 0.87636],
     [-2.601576, -4.030372, 1.313858],
     [-2.101008, -4.64449, 0.10382],
     [-4.098304, 3.811641, -1.565514],
     [-4.077992, 4.633016, -2.088623],
     [-4.820794, 4.037279, -0.937515],
     [-3.616287, 0.549545, -2.677343],
     [-3.564851, 0.0539, -1.87592],
     [-4.443122, 1.048319, -2.659152],
     [0.093079, 3.499639, -2.688856],
     [0.01739, 3.998067, -3.50861],
     [-0.841914, 3.248714, -2.53981],
     [-3.71763, 1.17797, -5.190931],
     [-3.468531, 0.398126, -5.73694],
     [-3.534477, 0.834868, -4.270742],
     [2.189021, -1.540005, 4.250915],
     [2.389576, -1.0943, 3.445708],
     [1.557192, -0.910271, 4.695886],
     [0.219511, -1.813941, -0.834816],
     [0.815111, -2.259853, -0.222439],
     [0.459479, -2.23041, -1.665345],
     [0.169965, -3.846934, 5.235632],
     [0.372012, -4.771048, 4.929771],
     [0.871485, -3.365644, 4.821651],
     [2.690342, -3.349566, 0.008587],
     [2.154892, -4.182669, 0.206891],
     [3.341973, -3.387489, 0.725831]]
)
FRAGMENTS = (
    [[0, 1, 2],
     [3, 4, 5],
     [6, 7, 8],
     [9, 10, 11],
     [12, 13, 14],
     [15, 16, 17],
     [18, 19, 20],
     [21, 22, 23],
     [24, 25, 26],
     [27, 28, 29],
     [30, 31, 32],
     [33, 34, 35],
     [36, 37, 38],
     [39, 40, 41],
     [42, 43, 44],
     [45, 46, 47]]
)
LEVEL_TO_SPEC = {2: 'mp2_aug_cc_pvdz_lvl2', 3: 'mp2_aug_cc_pvdz_lvl3'}
SPEC_KEYWORDS = {'d_convergence': 1e-06, 'e_convergence': 1e-06, 'maxiter': 100, 'scf_type': 'df'}
BSSE_TYPES = ['nocp']
MAX_NBODY = 4
RETURN_TOTAL_DATA = True
SUPERSYSTEM_IE_ONLY = False
TOTAL_CHARGE = 0
TOTAL_MULTIPLICITY = 1


def create_molecule() -> Molecule:
    """Construct the fragmented molecular system."""

    return Molecule(
        symbols=SYMBOLS,
        geometry=GEOMETRY,
        fragments=FRAGMENTS,
        molecular_charge=TOTAL_CHARGE,
        molecular_multiplicity=TOTAL_MULTIPLICITY,
    )


def create_specifications() -> Dict[str, AtomicSpecification]:
    """Build per-level atomic specifications."""

    return {
        label: AtomicSpecification(
            program="psi4",
            driver="energy",
            model={"method": "mp2", "basis": "aug-cc-pvdz"},
            keywords=dict(SPEC_KEYWORDS),
            protocols={"stdout": True},
            extras={},
        )
        for label in set(LEVEL_TO_SPEC.values())
    }


def create_manybody_input() -> ManyBodyInput:
    """Create the QCManyBody input model for the calculation."""

    molecule = create_molecule()
    keywords = ManyBodyKeywords(
        bsse_type="nocp",
        max_nbody=MAX_NBODY,
        return_total_data=RETURN_TOTAL_DATA,
        supersystem_ie_only=SUPERSYSTEM_IE_ONLY,
        levels=LEVEL_TO_SPEC,
    )

    specification = ManyBodySpecification(
        driver="energy",
        keywords=keywords,
        specification=create_specifications(),
    )

    return ManyBodyInput(molecule=molecule, specification=specification)


def create_parallel_config() -> ParallelConfig:
    """Return the parallel execution configuration."""

    return ParallelConfig(
        max_workers=6,
        execution_mode="multiprocessing",
        use_qcengine=True,
        qc_program="psi4",
        basis_set="aug-cc-pvdz",
        memory_limit_mb=1000,
        timeout_seconds=7200,
        qcengine_config={
            "keywords": dict(SPEC_KEYWORDS),
            "protocols": {"stdout": False},
        },
    )


def run_parallel_calculation():
    """Execute the many-body calculation using the parallel executor."""

    manybody_input = create_manybody_input()
    config = create_parallel_config()
    executor = ParallelManyBodyExecutor.from_manybodyinput(
        manybody_input,
        config,
    )

    start = time.time()
    fragment_results = executor.execute_full_calculation()
    stats = executor.get_execution_statistics()
    elapsed = time.time() - start

    print("QCManyBody parallel calculation complete!")
    total_fragments = stats.get("total_fragments")
    if total_fragments is not None:
        print(f"Total fragments executed: {total_fragments}")
    levels_executed = stats.get("levels_executed")
    if levels_executed is not None:
        print(f"Levels processed: {levels_executed}")
    parallel_time = stats.get("parallel_time")
    if parallel_time is not None:
        print(f"Parallel wall time: {parallel_time:.2f} s")
    serial_time = stats.get("serial_time")
    if serial_time:
        print(f"Serial baseline: {serial_time:.2f} s")
    speedup_factor = stats.get("speedup_factor")
    if speedup_factor:
        print(f"Measured speedup: {speedup_factor:.2f}x")
    print(f"Elapsed (wall clock): {elapsed:.2f} s")

    analysis = executor.core.analyze(fragment_results)
    results = analysis.get("results", {}) if isinstance(analysis, dict) else {}
    if results:
        total_key = f"nocp_corrected_total_energy_through_{MAX_NBODY}_body"
        interaction_key = (
            f"nocp_corrected_interaction_energy_through_{MAX_NBODY}_body"
        )
        total_energy = results.get(total_key)
        interaction_energy = results.get(interaction_key)
        if total_energy is not None:
            print(
                f"Total energy (through {MAX_NBODY}-body): "
                f"{float(total_energy):.12f} Eh"
            )
        if interaction_energy is not None:
            print(
                f"Interaction energy (through {MAX_NBODY}-body): "
                f"{float(interaction_energy):.12f} Eh"
            )

    return analysis


if __name__ == "__main__":
    run_parallel_calculation()
