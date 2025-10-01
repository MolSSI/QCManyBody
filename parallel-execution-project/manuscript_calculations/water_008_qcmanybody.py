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
    [[-1.266782, 0.1573, 2.395284],
     [-1.460012, -0.53144, 1.739857],
     [-0.989983, -0.313128, 3.206585],
     [-3.982693, -3.592347, 1.877175],
     [-4.513793, -2.843146, 1.957387],
     [-4.201233, -4.125344, 2.65745],
     [1.644004, -1.309028, 4.562366],
     [2.007636, -1.292204, 5.431114],
     [2.347273, -0.851564, 4.023934],
     [-3.339777, 1.183483, 4.633492],
     [-2.606423, 1.587589, 5.150428],
     [-3.978816, 0.94231, 5.362517],
     [0.181272, -3.550094, -4.962615],
     [0.139425, -4.178544, -5.731809],
     [0.019807, -2.72576, -5.397992],
     [-2.449372, -1.103753, -0.35684],
     [-2.684012, -2.006143, -0.348781],
     [-1.599227, -1.139902, -0.688068],
     [1.283864, -3.089752, 1.033513],
     [0.58277, -3.683348, 1.453278],
     [1.596638, -2.610531, 1.816469],
     [-0.030328, -1.919092, -1.36276],
     [0.138918, -2.478884, -0.59688],
     [0.671576, -2.196301, -1.955461],
     [0.020594, 3.044418, -2.922689],
     [0.163998, 3.923658, -3.286732],
     [0.756462, 2.550836, -3.340137],
     [0.100994, 2.880509, 3.083538],
     [-0.411779, 2.98813, 2.265206],
     [0.060619, 1.931607, 3.21503],
     [-3.364471, 1.078093, -2.794622],
     [-3.272216, 0.856308, -1.881994],
     [-2.593696, 1.602681, -3.046563],
     [-0.967763, 3.737534, 0.528569],
     [-1.888037, 3.494379, 0.48782],
     [-0.889402, 4.078984, -0.379858],
     [3.74307, -0.987074, 1.418783],
     [3.556224, -0.017249, 1.391744],
     [4.545072, -1.017973, 0.957457],
     [3.852465, 4.061565, -0.84645],
     [2.951219, 4.420497, -0.933834],
     [4.158999, 4.580643, -0.069352],
     [3.858351, 3.503932, 4.27209],
     [4.302934, 2.663119, 4.58],
     [3.213021, 3.186333, 3.581908],
     [3.421712, 1.478166, 0.56498],
     [2.835682, 2.245826, 0.733941],
     [3.890247, 1.874254, -0.178004]]
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
        max_workers=12,
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
        for n_body in range(1, MAX_NBODY + 1):
            n_body_key = f"nocp_{n_body}_body_energy"
            n_body_energy = results.get(n_body_key)
            if n_body_energy is not None:
                print(f"{n_body}-body energy: {float(n_body_energy):.12f} Eh")
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
