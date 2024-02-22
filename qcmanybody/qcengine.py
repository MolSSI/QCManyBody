from typing import Mapping, Literal, Union, Iterable, Any
from pprint import pprint

import qcengine as qcng
from qcelemental.models import AtomicInput

from qcmanybody.manybody import Molecule, ManyBodyCalculator
from qcmanybody.models import BsseEnum
from qcmanybody.utils import delabeler


def run_qcengine(
        molecule: Molecule,
        levels: Mapping[Union[int, Literal["supersystem"]], str],
        specifications: Mapping[str, Mapping[str, Any]],
        bsse_type: Iterable[BsseEnum],
        return_total_data: bool,
):

    mc = ManyBodyCalculator(molecule, bsse_type, levels, return_total_data)

    component_results = {}

    computation_count = {}
    for chem, label, imol in mc.iterate_molecules(True):
        print(label)
        inp = AtomicInput(molecule=imol, **specifications[chem]['specification'])

        _, real, bas = delabeler(label)
        computation_count.setdefault(len(real), 0)
        computation_count[len(real)] += 1

        result = qcng.compute(inp, specifications[chem]["program"])

        if not result.success:
            print(result.error.error_message)
            raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

        component_results[label] = {
            'energy': result.properties.return_energy,
            'gradient': result.properties.return_gradient,
            'hessian': result.properties.return_hessian
        }


    print("COMPUTATION COUNTS")
    pprint(computation_count)


    return mc.analyze(component_results)
