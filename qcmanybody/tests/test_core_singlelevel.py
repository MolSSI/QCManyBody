import pytest
import qcelemental

from qcmanybody import ManyBodyCore
from qcmanybody.models import BsseEnum

from .common import mol_h2o_3_dict
from .utils import compare_results, load_component_data, load_ref_data


@pytest.mark.parametrize(
    "levels, component_file, ref_file",
    [
        # fmt: off
        # Grad reference file includes energies
        ({1: "e_scf", 2: "e_scf", 3: "e_scf"}, "h2o_trimer_single_energy_1",   "h2o_trimer_single_gradient_1"),
        ({1: "e_scf", 2: "e_scf",           }, "h2o_trimer_single_energy_2",   "h2o_trimer_single_gradient_2"),
        ({1: "g_scf", 2: "g_scf", 3: "g_scf"}, "h2o_trimer_single_gradient_1", "h2o_trimer_single_gradient_1"),
        ({1: "g_scf", 2: "g_scf"            }, "h2o_trimer_single_gradient_2", "h2o_trimer_single_gradient_2"),
        ({1: "h_scf", 2: "h_scf", 3: "h_scf"}, "h2o_trimer_single_hessian_1",  "h2o_trimer_single_hessian_1"),
        ({1: "h_scf", 2: "h_scf"            }, "h2o_trimer_single_hessian_2",  "h2o_trimer_single_hessian_2"),
        # fmt: on
    ],
)
def test_h2o_trimer_single(levels, component_file, ref_file):

    component_results = load_component_data(component_file)
    ref_data = load_ref_data(ref_file)

    mc = ManyBodyCore(mol_h2o_3_dict, [BsseEnum.cp, BsseEnum.nocp, BsseEnum.vmfc], levels,
        embedding_charges=None, supersystem_ie_only=False, return_total_data=True)
    nbody_results = mc.analyze(component_results)
    compare_results(nbody_results, ref_data, levels)


def test_core_mol_error():
    # check sensible error from internal Molecule construction
    odd_mol = mol_h2o_3_dict.copy()
    odd_mol["symbols"][0] = "A"
    with pytest.raises(qcelemental.exceptions.NotAnElementError):
        ManyBodyCore(odd_mol, [BsseEnum.nocp], {1: "asdf"},
            return_total_data=False, supersystem_ie_only=False, embedding_charges=None)
