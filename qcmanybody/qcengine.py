from typing import Mapping, Literal, Union, Iterable, Any

import qcengine as qcng
from qcelemental.models import AtomicInput

from qcmanybody.manybody import Molecule, ManybodyCalculator
from qcmanybody.models import BsseEnum


def run_qcengine(
        molecule: Molecule,
        levels: Mapping[Union[int, Literal["supersystem"]], str],
        specifications: Mapping[str, Mapping[str, Any]],
        bsse_type: Iterable[BsseEnum],
        return_total_data: bool,
):

    mc = ManybodyCalculator(molecule, bsse_type, levels, return_total_data)

    component_results = {}

    for chem, label, imol in mc.iterate_molecules(True):
        inp = AtomicInput(molecule=imol, **specifications[chem]['specification'])
        result = qcng.compute(inp, specifications[chem]["program"])

        if not result.success:
            print(result.error.error_message)
            raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

        component_results[label] = {
            'energy': result.properties.return_energy,
            'gradient': result.properties.return_gradient,
            'hessian': result.properties.return_hessian
        }

    print(component_results)
    return mc.analyze(component_results)