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
    [[-1.314959, -0.355371, 1.609453],
     [-1.650552, -0.466153, 0.705909],
     [-1.904209, 0.297125, 2.037428],
     [-2.576046, 2.52114, 4.576931],
     [-2.03952, 1.730785, 4.813177],
     [-1.883294, 3.240249, 4.608381],
     [-3.823533, 2.967636, -0.174021],
     [-3.071116, 3.110859, 0.438056],
     [-4.409235, 2.572923, 0.481411],
     [0.456002, 2.799863, 3.55061],
     [0.862763, 1.956281, 3.291571],
     [0.376297, 3.245106, 2.705187],
     [1.451225, -2.844857, 0.731514],
     [1.175653, -3.711989, 0.383871],
     [1.432135, -2.322698, -0.101716],
     [-3.840783, -3.293364, 2.407645],
     [-4.302789, -2.503474, 2.296166],
     [-4.482447, -3.889175, 2.824719],
     [-1.1264, -3.895366, 1.886693],
     [-2.099847, -3.645839, 1.987824],
     [-0.708043, -3.151914, 2.347949],
     [0.845414, -2.100188, 3.072418],
     [1.360598, -2.492743, 2.388517],
     [0.598455, -1.218188, 2.678897],
     [-3.246965, -1.514128, -0.514774],
     [-2.61138, -2.186002, -0.396261],
     [-4.014363, -2.002173, -0.596489],
     [1.133182, -3.619761, -4.971056],
     [1.217079, -4.287052, -5.703208],
     [0.682214, -2.919569, -5.419953],
     [-2.762656, 1.116332, -2.249623],
     [-3.168987, 1.637842, -1.576191],
     [-2.203618, 1.702895, -2.775154],
     [-0.496519, -1.957706, -1.535191],
     [-0.219896, -2.872669, -1.657279],
     [-0.280062, -1.581385, -2.390958],
     [3.962929, 2.405329, 2.481884],
     [4.655379, 2.506414, 3.195833],
     [3.482595, 1.570992, 2.740404],
     [1.222824, 1.38071, 0.432964],
     [0.812804, 0.661211, 0.904019],
     [2.01378, 1.462024, 0.994878],
     [0.881979, 3.882261, -1.312077],
     [0.718444, 3.008846, -0.942541],
     [1.815579, 3.790349, -1.593735],
     [3.476071, -0.351383, 1.265784],
     [2.958736, 0.230804, 0.657816],
     [4.155388, -0.626003, 0.700013]]
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
