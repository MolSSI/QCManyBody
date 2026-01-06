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
    [[3.366823, 1.954801, 2.137534],
     [3.979546, 2.541644, 2.66632],
     [3.355147, 1.103634, 2.656233],
     [1.996464, 1.260638, -0.371882],
     [2.493682, 1.800178, 0.235869],
     [2.667157, 0.566306, -0.498548],
     [3.596148, -1.148118, 1.764902],
     [2.893943, -1.462762, 1.145136],
     [4.354434, -1.455145, 1.331634],
     [1.500647, 4.146142, -1.141178],
     [1.313737, 3.231897, -0.905857],
     [2.428898, 4.227127, -0.839215],
     [1.967982, -2.051116, -0.078714],
     [1.979806, -1.538294, -0.906716],
     [1.251951, -2.690046, -0.294015],
     [-0.092643, -2.462756, -2.572997],
     [0.057498, -2.965036, -3.38155],
     [0.26595, -3.059844, -1.912918],
     [4.160926, -4.228708, -0.706753],
     [4.276987, -4.437521, 0.244146],
     [4.890941, -3.602002, -0.759242],
     [0.874813, -3.901487, -4.487329],
     [1.671277, -4.115955, -3.932341],
     [1.044123, -4.431107, -5.252827],
     [-3.765549, -1.355387, -1.48374],
     [-3.141744, -1.611961, -0.83995],
     [-3.448894, -0.534091, -1.726572],
     [-2.049267, 3.309727, 1.619404],
     [-2.616326, 4.092907, 1.911235],
     [-1.549315, 3.144675, 2.433839],
     [-2.043356, 2.974812, -4.638324],
     [-2.316972, 2.0384, -4.509488],
     [-1.143148, 2.856534, -5.055036],
     [-3.328697, 1.674209, -1.26345],
     [-3.091867, 2.52802, -1.588269],
     [-3.752934, 1.796391, -0.404464],
     [-0.376536, -0.681051, 0.970989],
     [0.2468, -0.048589, 0.580182],
     [-1.221373, -0.555613, 0.494917],
     [-4.122445, -3.788983, 2.815779],
     [-4.321143, -2.908578, 3.003418],
     [-4.772761, -4.295353, 3.327006],
     [0.402133, 2.128508, 4.26094],
     [0.123101, 1.896623, 3.359521],
     [1.001898, 1.41017, 4.469766],
     [0.484762, -0.608337, 3.705584],
     [0.275069, -1.379411, 3.206878],
     [-0.235654, -0.592957, 4.394449]]
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
