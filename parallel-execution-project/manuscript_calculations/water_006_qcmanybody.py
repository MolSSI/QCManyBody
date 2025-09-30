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
    [[3.68185, 2.794538, 2.87362],
     [4.384392, 2.652901, 3.570615],
     [3.165972, 1.941682, 2.886895],
     [-0.053172, 3.157496, 3.266603],
     [0.717892, 3.414325, 2.733982],
     [-0.635022, 2.787117, 2.600566],
     [3.610915, -0.105882, 1.14047],
     [3.524493, 0.806244, 0.770659],
     [4.436671, -0.349155, 0.799986],
     [0.590179, 1.882994, 0.637364],
     [0.315494, 1.757073, -0.266175],
     [0.618286, 0.943263, 0.890515],
     [0.837336, -3.628563, -4.785764],
     [0.320134, -4.419663, -5.094024],
     [0.941009, -3.156434, -5.599088],
     [1.921641, -3.324268, 0.489865],
     [1.180185, -3.856223, 0.830424],
     [1.419522, -2.561692, 0.124298],
     [-0.192274, -2.240961, -1.377005],
     [-0.300864, -3.125943, -1.011494],
     [-0.143242, -2.431526, -2.316201],
     [-0.833612, -3.149459, 1.321281],
     [-1.643781, -3.488354, 1.820147],
     [-0.131056, -3.399989, 1.941079],
     [-1.038049, -0.280804, 1.997115],
     [-0.382953, -0.880978, 2.386893],
     [-1.77891, -0.846288, 1.701345],
     [-2.425465, 2.065041, 4.506157],
     [-2.801669, 1.675718, 5.32788],
     [-1.480343, 1.745413, 4.557121],
     [-3.414001, -3.584178, 2.350093],
     [-3.691525, -2.84867, 1.868624],
     [-4.21641, -3.893299, 2.798711],
     [0.73216, -1.67768, 3.534502],
     [1.475958, -1.422415, 3.016006],
     [0.409495, -0.806013, 3.894856],
     [-3.830254, 2.057838, -0.157201],
     [-3.141579, 1.769266, 0.478199],
     [-4.132481, 2.811528, 0.361495],
     [0.195532, 3.55846, -2.177552],
     [0.538395, 3.13683, -1.383304],
     [1.027023, 3.739625, -2.662518],
     [-2.849135, -1.48901, -0.018304],
     [-3.439426, -1.959073, -0.56605],
     [-2.110085, -2.024567, -0.045657],
     [-2.693852, 1.234169, -2.488064],
     [-3.317738, 1.520489, -1.840472],
     [-1.828107, 1.575626, -2.229886]]
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
