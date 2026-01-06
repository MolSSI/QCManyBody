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
    [[2.005979, -4.685521, -1.844031],
     [1.355527, -4.283326, -1.268374],
     [1.526601, -4.878857, -2.649665],
     [2.608892, -2.675449, 2.890175],
     [2.427245, -3.533869, 2.507619],
     [2.493361, -2.060559, 2.16575],
     [0.069765, -3.765074, -0.092294],
     [-0.748399, -4.102246, -0.457193],
     [-0.205709, -3.169585, 0.604654],
     [0.0, 0.0, 0.0],
     [-0.900139, -0.263684, 0.190917],
     [0.092177, -0.127148, -0.944228],
     [2.017875, -5.101032, 1.831973],
     [2.301866, -5.927774, 1.442003],
     [1.219745, -4.872587, 1.355497],
     [-2.841906, -3.401469, -1.848014],
     [-2.92869, -2.964138, -1.000996],
     [-2.511626, -2.722254, -2.43607],
     [3.585351, -2.455495, -1.561401],
     [2.94275, -3.137891, -1.75538],
     [4.413623, -2.925805, -1.466523],
     [2.450667, -1.324316, 0.643786],
     [2.758818, -1.643552, -0.204364],
     [1.562815, -1.009682, 0.473627],
     [1.078891, 0.416333, 5.43754],
     [1.696153, -0.313147, 5.493008],
     [1.275802, 0.834104, 4.599135],
     [-2.796958, -0.664606, -0.330683],
     [-3.203523, 0.066005, 0.135301],
     [-2.926777, -0.460567, -1.256827],
     [-3.815253, -2.427004, 2.120311],
     [-4.683296, -2.813738, 2.00556],
     [-3.776899, -1.726033, 1.46962],
     [-1.517278, 1.84158, 2.370512],
     [-2.091013, 1.341235, 2.950781],
     [-1.58883, 1.400826, 1.523844],
     [-1.047178, -2.933996, 2.272746],
     [-0.552524, -2.353539, 2.851206],
     [-1.964328, -2.754123, 2.479411],
     [-4.656323, 1.08899, 0.825096],
     [-4.978444, 1.694752, 0.157623],
     [-5.185808, 1.281009, 1.599048],
     [-3.517569, -1.259768, 4.504909],
     [-4.256747, -0.824819, 4.929953],
     [-3.907261, -1.729748, 3.767692],
     [-0.799077, -1.047125, 4.089732],
     [-0.345323, -0.468214, 4.702266],
     [-1.72862, -0.879239, 4.244624],
     [4.160764, 2.166572, 2.715175],
     [4.287608, 3.015355, 2.291257],
     [3.210581, 2.05839, 2.756142],
     [-1.138153, 4.410088, 0.263741],
     [-0.296342, 4.68125, 0.629877],
     [-1.763663, 4.543124, 0.975969],
     [3.185867, 3.687994, -0.593526],
     [3.28441, 4.021231, -1.485418],
     [2.23966, 3.656005, -0.452463],
     [-0.359419, 3.816258, 3.849684],
     [0.347573, 4.414256, 3.60722],
     [-0.624062, 3.410303, 3.024218],
     [-3.090156, 3.671481, -1.577749],
     [-2.493968, 3.850117, -0.850509],
     [-3.963404, 3.713109, -1.187962],
     [1.502374, 1.25824, 2.86289],
     [1.112074, 0.55707, 2.34111],
     [1.041895, 2.049705, 2.584017],
     [0.604867, 4.681388, -1.883436],
     [0.594428, 5.575712, -2.224471],
     [-0.212207, 4.60275, -1.39106],
     [3.780276, 0.979797, -0.108064],
     [3.615404, 0.632822, 0.768666],
     [3.949462, 1.911673, 0.030538],
     [-2.007336, 1.839099, -5.027128],
     [-2.54765, 2.629216, -5.024783],
     [-1.19034, 2.105807, -4.605689],
     [3.379245, 2.098633, -4.312302],
     [3.524446, 2.9781, -3.963471],
     [3.114449, 1.579866, -3.552701],
     [-2.060018, -1.466237, -3.555814],
     [-1.262512, -0.993527, -3.317574],
     [-1.848445, -1.889886, -4.387672],
     [1.88149, 0.12892, -6.276787],
     [2.830236, 0.06493, -6.167178],
     [1.537069, -0.655066, -5.849038],
     [-3.987328, 0.453972, -2.785206],
     [-3.613701, 1.3271, -2.904706],
     [-3.470783, -0.10907, -3.361742],
     [0.301391, -0.222093, -2.917628],
     [1.235922, -0.426827, -2.886593],
     [0.264685, 0.687556, -3.213298],
     [3.043168, 0.112457, -2.60588],
     [3.945524, -0.17369, -2.464093],
     [2.730091, 0.359177, -1.735628],
     [0.568276, 2.339609, -3.96175],
     [1.254414, 2.578278, -4.585032],
     [0.66538, 2.967531, -3.24585]]
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
     [45, 46, 47],
     [48, 49, 50],
     [51, 52, 53],
     [54, 55, 56],
     [57, 58, 59],
     [60, 61, 62],
     [63, 64, 65],
     [66, 67, 68],
     [69, 70, 71],
     [72, 73, 74],
     [75, 76, 77],
     [78, 79, 80],
     [81, 82, 83],
     [84, 85, 86],
     [87, 88, 89],
     [90, 91, 92],
     [93, 94, 95]]
)
LEVEL_TO_SPEC = {
    1: "mp2_aug_cc_pvdz_lvl1",
    2: "mp2_aug_cc_pvdz_lvl2",
    3: "mp2_aug_cc_pvdz_lvl3",
    4: "mp2_aug_cc_pvdz_lvl4",
}
SPEC_KEYWORDS = {
    "d_convergence": 1e-06,
    "e_convergence": 1e-06,
    "maxiter": 100,
    "scf_type": "df",
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
