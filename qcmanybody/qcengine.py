from typing import Mapping, Literal, Union, Iterable, Any

import qcengine as qcng
from qcelemental.models import AtomicInput

from qcmanybody.manybody import Molecule, ManyBodyCalculator
from qcmanybody.models import BsseEnum
from qcmanybody.utils import delabeler


def run_qcengine_base(
    molecule: Molecule,
    levels: Mapping[Union[int, Literal["supersystem"]], str],
    specifications: Mapping[str, Mapping[str, Any]],
    bsse_type: Iterable[BsseEnum],
    return_total_data: bool,
):

    mc = ManyBodyCalculator(molecule, bsse_type, levels, return_total_data)

    component_results = {}

    computation_count = {}
    for chem, label, imol in mc.iterate_molecules():
        print(label)
        inp = AtomicInput(molecule=imol, **specifications[chem]["specification"])

        _, real, bas = delabeler(label)
        computation_count.setdefault(len(real), 0)
        computation_count[len(real)] += 1

        result = qcng.compute(inp, specifications[chem]["program"])

        if not result.success:
            print(result.error.error_message)
            raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

        # pull out stuff
        props = {"energy", "gradient", "hessian"}

        component_results[label] = {}

        for p in props:
            if hasattr(result.properties, f"return_{p}"):
                v = getattr(result.properties, f"return_{p}")
                # print(f"  {label} {p}: {v}")
                if v is not None:
                    component_results[label][p] = v

    return mc, component_results


def run_qcengine(
    molecule: Molecule,
    levels: Mapping[Union[int, Literal["supersystem"]], str],
    specifications: Mapping[str, Mapping[str, Any]],
    bsse_type: Iterable[BsseEnum],
    return_total_data: bool,
):

    mc, component_results = run_qcengine_base(molecule, levels, specifications, bsse_type, return_total_data)
    return mc.analyze(component_results)
