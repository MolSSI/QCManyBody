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
    [[-1.495337, 0.711664, 2.282648],
     [-1.089828, 0.433028, 1.44646],
     [-2.422223, 0.402049, 2.248732],
     [0.025758, 3.464002, 2.206748],
     [-0.474684, 4.134329, 1.712376],
     [-0.41909, 2.656237, 1.944144],
     [-2.394802, 4.25259, 0.743078],
     [-3.074574, 4.656795, 1.27432],
     [-2.9255, 3.508241, 0.408004],
     [-3.998369, -0.409021, 2.000685],
     [-4.894572, -0.473348, 1.589731],
     [-3.492098, -0.867122, 1.37553],
     [0.055698, -4.108929, -2.434314],
     [0.09709, -4.310531, -3.37442],
     [-0.876752, -4.336368, -2.238872],
     [3.108913, -3.572295, -0.612321],
     [2.326801, -4.21124, -0.600429],
     [2.998494, -3.147014, 0.252227],
     [3.701016, -3.343872, 2.143097],
     [3.829909, -2.503357, 2.499081],
     [3.348073, -3.864862, 2.881118],
     [0.69652, -1.899657, -0.663703],
     [1.363919, -1.737683, -1.339668],
     [0.517655, -2.832858, -0.797682],
     [-2.114587, -1.381906, -0.733854],
     [-1.886865, -2.129508, -0.225264],
     [-1.305396, -0.967916, -0.820951],
     [-3.872933, 1.117779, -4.931776],
     [-3.111494, 1.557278, -4.489754],
     [-4.231177, 0.557832, -4.186032],
     [-3.165684, 3.60364, -2.461012],
     [-3.982203, 3.075836, -2.519593],
     [-3.429211, 4.231644, -1.751469],
     [-3.929118, 0.484813, -2.282651],
     [-4.359475, -0.170564, -1.757424],
     [-2.987632, 0.457628, -2.069043],
     [3.12936, 3.42452, 4.416857],
     [3.797435, 2.889631, 4.933576],
     [2.323395, 2.838112, 4.402159],
     [1.609917, 0.909086, 2.406862],
     [1.307928, 1.840765, 2.452238],
     [0.728223, 0.520472, 2.401364],
     [2.525879, -1.690961, 3.936806],
     [2.229646, -0.803417, 3.828431],
     [3.350173, -1.578985, 4.486144],
     [0.640689, 4.156085, 5.056224],
     [0.649532, 3.763687, 4.142819],
     [0.900493, 5.046444, 4.869372]]
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
