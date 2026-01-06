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
    [[1.461, 2.754, -3.057],
     [2.412, 2.627, -2.982],
     [1.209, 1.956, -3.566],
     [-1.722, -0.078, -3.8],
     [-1.299, -0.608, -4.44],
     [-2.444, -0.594, -3.585],
     [-1.579, -0.225, 0.058],
     [-1.316, -1.066, -0.348],
     [-0.781, 0.34, 0.048],
     [-2.339, 2.681, -3.057],
     [-2.057, 1.739, -3.019],
     [-3.222, 2.599, -3.517],
     [2.14, -0.331, -1.296],
     [2.732, -0.163, -2.051],
     [2.807, -0.592, -0.622],
     [-0.68, -2.561, -1.307],
     [-0.982, -2.737, -2.205],
     [0.218, -2.263, -1.467],
     [3.572, -3.58, -1.522],
     [3.636, -3.661, -2.497],
     [2.628, -3.766, -1.47],
     [1.141, -2.389, -4.501],
     [0.697, -2.89, -3.766],
     [2.031, -2.698, -4.414],
     [-3.551, 3.366, 4.14],
     [-3.798, 2.479, 4.095],
     [-3.204, 3.563, 3.256],
     [-1.29, 2.185, 2.54],
     [-1.533, 2.236, 1.561],
     [-0.59, 2.855, 2.58],
     [-3.312, -1.219, 2.114],
     [-2.512, -0.81, 1.703],
     [-3.947, -0.571, 1.93],
     [-2.867, 3.08, 0.078],
     [-2.799, 3.102, -0.863],
     [-3.642, 2.549, 0.302],
     [0.654, 1.717, 0.023],
     [1.059, 1.126, -0.605],
     [1.409, 1.785, 0.634],
     [0.539, -1.395, 3.295],
     [1.46, -1.298, 3.123],
     [0.434, -2.384, 3.363],
     [2.439, 1.144, 2.455],
     [3.31, 0.836, 2.837],
     [2.121, 1.819, 3.116],
     [-0.597, 0.687, 4.815],
     [-0.891, 1.272, 4.097],
     [-0.183, -0.027, 4.327]]
)

SYMBOLS = SYMBOLS[: len(GEOMETRY)]
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
LEVEL_TO_SPEC = {
    1: 'mp2_aug_cc_pvdz_lvl1',
    2: 'mp2_aug_cc_pvdz_lvl2',
    3: 'mp2_aug_cc_pvdz_lvl3',
    4: 'mp2_aug_cc_pvdz_lvl4',
}
SPEC_KEYWORDS = {
    'd_convergence': 1e-06,
    'e_convergence': 1e-06,
    'maxiter': 100,
    'scf_type': 'df',
}
BSSE_TYPES = "nocp"
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
            n_body_key = f"nocp_corrected_total_energy_through_{n_body}_body"
            n_body_energy = results.get(n_body_key)
            if n_body_energy is not None:
                print(f"{n_body}-body energy: {float(n_body_energy):.12f} Eh")
                contribution_key = (
                    "nocp_corrected_"
                    f"{n_body}_body_contribution_to_energy"
                )
                contribution = results.get(contribution_key)
                if contribution is not None:
                    print(
                        "  "
                        f"Contribution from {n_body}-body order: "
                        f"{float(contribution):.12f} Eh"
                    )
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
