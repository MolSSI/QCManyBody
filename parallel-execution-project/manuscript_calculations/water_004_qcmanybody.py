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
    [[1.289545, 3.622126, -1.234222],
     [1.516741, 2.736353, -0.934317],
     [2.006483, 3.786236, -1.881103],
     [-0.623906, -0.327068, 1.209419],
     [-0.552412, 0.228582, 0.417315],
     [-1.5621, -0.59689, 1.265256],
     [-2.960615, 1.551458, -1.693306],
     [-2.407707, 2.064065, -2.26086],
     [-3.33713, 2.149564, -1.035076],
     [1.294444, 1.23841, 0.012121],
     [1.076204, 1.09, 0.927564],
     [1.825869, 0.433479, -0.120713],
     [-2.122137, 3.447414, 1.892218],
     [-2.826831, 4.051366, 2.29064],
     [-1.317729, 3.890921, 2.20323],
     [-2.225668, 2.814883, 4.904774],
     [-2.324615, 1.967662, 5.395445],
     [-1.415078, 3.198525, 5.344975],
     [0.656126, 2.601905, 4.050683],
     [0.830957, 1.865698, 4.660287],
     [-0.119244, 2.287033, 3.582716],
     [-3.779167, 3.893046, -0.392587],
     [-2.954598, 4.248679, 0.000977],
     [-4.321495, 4.008296, 0.39547],
     [1.517185, -2.080446, 0.275877],
     [1.370479, -2.993096, -0.031144],
     [0.767312, -1.634289, -0.177908],
     [1.157642, -3.919991, -4.712879],
     [1.482782, -4.597069, -5.364193],
     [1.125355, -3.148726, -5.259923],
     [-3.374181, -1.753947, -1.126278],
     [-3.632256, -1.032976, -1.658281],
     [-3.378912, -2.437496, -1.731673],
     [-0.169102, -2.277746, -1.880495],
     [-0.59266, -3.129658, -2.033563],
     [-0.876898, -1.672861, -2.112759],
     [3.346885, -0.739109, 1.455676],
     [2.495949, -1.146635, 1.162364],
     [3.954914, -1.284332, 1.019766],
     [0.635713, -1.390592, 3.660639],
     [1.546539, -1.624232, 3.715803],
     [0.562533, -1.01904, 2.738487],
     [4.023571, 2.373326, 2.358002],
     [4.980223, 2.433122, 2.642029],
     [3.944691, 1.446061, 2.000753],
     [4.15173, -3.364424, 2.461883],
     [3.843025, -2.655099, 2.963243],
     [3.347114, -3.828962, 2.183542]]
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
