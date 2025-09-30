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
    [[-4.167207, 1.760194, 0.11044],
     [-3.409954, 2.327637, 0.36707],
     [-4.789634, 2.483055, -0.025545],
     [-0.26248, 2.835114, -2.840283],
     [0.190162, 2.918128, -3.685492],
     [-1.084127, 3.336273, -3.022283],
     [-0.235445, 2.547268, 0.657216],
     [0.112008, 2.25486, -0.180319],
     [-0.142717, 1.702179, 1.131764],
     [-2.576816, 1.575475, -2.962402],
     [-2.714304, 1.43154, -2.039922],
     [-1.624698, 1.647275, -3.107582],
     [0.300694, -3.819203, -4.933036],
     [0.159399, -4.246794, -5.819357],
     [0.875648, -3.106064, -5.169704],
     [-0.145604, -2.980895, 1.317655],
     [-0.785926, -3.758038, 1.239389],
     [0.572947, -3.398015, 1.81785],
     [2.809315, -3.803742, -0.244227],
     [2.050744, -3.328262, 0.139442],
     [3.468817, -3.660921, 0.471277],
     [-0.006396, -2.114211, -1.52385],
     [0.643706, -1.908126, -0.843057],
     [0.397966, -2.876537, -1.943562],
     [-3.598599, -3.368084, 2.023607],
     [-3.334535, -2.659739, 2.551167],
     [-4.332037, -3.771225, 2.513731],
     [-1.094601, -0.048243, 2.084545],
     [-1.887154, 0.128996, 2.615329],
     [-1.345363, -0.751045, 1.452631],
     [-2.892897, 1.55136, 4.750626],
     [-3.203435, 1.01273, 5.513368],
     [-2.263246, 2.181002, 5.203509],
     [-2.404242, -1.294073, -0.069406],
     [-2.899232, -1.949582, -0.51068],
     [-1.563649, -1.456493, -0.386877],
     [3.687372, 3.071726, 3.542292],
     [4.455071, 2.43183, 3.517774],
     [3.419507, 3.076719, 4.502445],
     [3.256831, -0.801148, 1.404839],
     [3.214807, 0.170418, 1.230213],
     [4.129479, -0.980825, 1.153471],
     [0.112397, 3.299059, 3.351725],
     [-0.634815, 3.308272, 2.730617],
     [0.407322, 2.388969, 3.287752],
     [1.37509, -1.288864, 4.423921],
     [0.528697, -0.999088, 4.129162],
     [1.995922, -0.793278, 3.821648]]
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
