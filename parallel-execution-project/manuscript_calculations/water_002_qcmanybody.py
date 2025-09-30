#!/usr/bin/env python3
from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Dict

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from qcelemental.models import Molecule  # noqa: E402
from qcmanybody.models import (  # noqa: E402
    AtomicSpecification,
    ManyBodyInput,
    ManyBodyKeywords,
    ManyBodySpecification,
)
from qcmanybody.parallel import (  # noqa: E402
    ParallelConfig,
    ParallelManyBodyExecutor,
)

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
    [[-3.907819, -4.170752, 3.453443],
     [-4.322974, -3.501216, 2.974724],
     [-3.126272, -3.742276, 3.83589],
     [-3.479767, -1.338075, -1.687071],
     [-3.42577, -1.031531, -0.808126],
     [-2.811748, -0.861264, -2.087279],
     [-0.644461, -0.489171, 0.653661],
     [-0.88867, 0.449128, 0.618408],
     [-1.173338, -0.931285, -0.039847],
     [-3.866717, -1.39377, 2.068426],
     [-3.727784, -0.448569, 1.816448],
     [-4.780618, -1.470139, 1.942182],
     [2.219004, -1.335169, -0.735943],
     [2.053367, -1.788545, -1.581949],
     [1.302267, -1.315756, -0.380273],
     [0.958356, -3.805365, -4.854227],
     [1.622768, -4.21664, -4.239596],
     [0.856572, -4.495707, -5.493146],
     [3.965362, -3.737, -1.315192],
     [4.460306, -3.206387, -1.974558],
     [3.647574, -2.985211, -0.803078],
     [-0.043444, -2.306831, -2.053961],
     [-0.179345, -2.145925, -2.994293],
     [-0.540335, -3.118277, -1.929615],
     [0.711654, 3.310394, -2.027755],
     [1.056294, 2.552213, -1.545532],
     [1.351114, 4.003017, -1.76176],
     [-3.358026, 1.924561, -0.756587],
     [-4.283849, 2.107429, -0.757274],
     [-3.070842, 1.87033, -1.677101],
     [-1.738787, 3.330576, 2.041959],
     [-2.678769, 3.57758, 2.316758],
     [-1.911847, 2.515912, 1.545096],
     [-1.705788, 2.877792, -3.809437],
     [-0.798685, 3.256433, -3.855559],
     [-1.906093, 2.739839, -4.778385],
     [1.960872, 1.354203, -0.392858],
     [2.323363, 0.492113, -0.574767],
     [1.824699, 1.232171, 0.563455],
     [0.426847, -0.839839, 3.447016],
     [0.408094, -1.418426, 2.703969],
     [0.736213, -1.4393, 4.180984],
     [0.103683, 2.099996, 4.411801],
     [-0.011235, 1.61912, 3.575298],
     [-0.282602, 1.489227, 5.041945],
     [2.779234, 1.603944, 2.293049],
     [3.130732, 2.524515, 2.461647],
     [2.218677, 1.418461, 3.096195]]
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
SPEC_KEYWORDS = {
    'd_convergence': 1e-06,
    'e_convergence': 1e-06,
    'maxiter': 100,
    'scf_type': 'df',
}
BSSE_TYPES = ['nocp']
MAX_NBODY = 3
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
