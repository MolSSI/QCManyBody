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
    [[-3.689394, 3.749647, 3.849786],
     [-3.619138, 4.643865, 3.637087],
     [-4.28168, 3.732044, 4.617615],
     [-2.81842, 2.504202, -0.275951],
     [-2.553654, 2.338385, -1.166452],
     [-3.558348, 1.916102, -0.077439],
     [-3.580505, -1.136058, 1.984962],
     [-2.658258, -1.05788, 1.639201],
     [-3.831792, -0.245211, 2.000118],
     [-1.465891, 2.685972, 2.15668],
     [-2.0156, 3.301605, 2.738837],
     [-1.859867, 2.86957, 1.289746],
     [1.190879, 2.896937, -2.871405],
     [0.942343, 2.639074, -1.978157],
     [1.913476, 3.531302, -2.684794],
     [-2.531796, -1.025968, -2.369001],
     [-2.151713, -1.064409, -3.219583],
     [-3.020945, -0.256827, -2.423033],
     [-2.162102, 2.865319, -2.935133],
     [-1.296293, 3.332618, -2.953543],
     [-2.139712, 2.38438, -3.810468],
     [-1.050886, -0.463277, 0.286565],
     [-1.387084, 0.421257, 0.500705],
     [-1.006015, -0.505299, -0.689319],
     [1.130504, 1.642189, -0.193717],
     [1.404608, 0.730152, -0.221036],
     [1.30475, 1.795748, 0.751813],
     [-0.32487, 1.518152, 4.266193],
     [0.050688, 1.801697, 3.41605],
     [-0.079652, 0.591471, 4.287789],
     [0.448535, -1.1353, 3.091499],
     [0.157273, -1.600171, 2.3258],
     [0.977792, -1.826755, 3.576831],
     [2.637379, 1.395002, 2.714204],
     [2.691847, 2.289005, 3.1583],
     [2.939143, 0.763914, 3.424371],
     [2.402853, -0.862671, -0.881986],
     [1.683715, -1.449266, -1.177724],
     [2.712186, -0.532529, -1.755259],
     [-0.305997, -2.575829, -1.775956],
     [0.196841, -2.171563, -1.0602],
     [-1.124501, -2.076707, -1.734276],
     [1.11661, -2.846958, -4.63215],
     [1.14411, -3.811086, -4.391191],
     [0.24109, -2.619637, -4.354769],
     [3.238682, -3.733328, -1.767495],
     [4.17876, -3.821178, -2.031726],
     [3.412336, -3.606945, -0.828182]]
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
