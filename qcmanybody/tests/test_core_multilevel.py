import pytest

from qcmanybody import ManyBodyCalculator  # test old name still operational
from qcmanybody.models import BsseEnum

from .common import mol_h2o_3
from .utils import compare_results, load_component_data, load_ref_data


@pytest.mark.parametrize(
    "levels, component_file, ref_file",
    [
        # fmt: off
        # Grad reference file includes energies
        ({1: "e_mp2", 2: "e_b3lyp", 3: "e_scf"}, "h2o_trimer_multi_energy_1",   "h2o_trimer_multi_gradient_1"),
        ({1: "e_mp2", 2: "e_b3lyp",           }, "h2o_trimer_multi_energy_2",   "h2o_trimer_multi_gradient_2"),
        ({1: "g_mp2", 2: "g_b3lyp", 3: "g_scf"}, "h2o_trimer_multi_gradient_1", "h2o_trimer_multi_gradient_1"),
        ({1: "g_mp2", 2: "g_b3lyp",           }, "h2o_trimer_multi_gradient_2", "h2o_trimer_multi_gradient_2"),
        # ({1: "h_scf_atz", 2: "h_scf_adz", 3: "h_scf"}, "h2o_trimer_multi_hessian_1", "h2o_trimer_multi_hessian_1"),
        # ({1: "h_scf_atz", 2: "h_scf_adz"            }, "h2o_trimer_multi_hessian_2", "h2o_trimer_multi_hessian_2"),
        # fmt: on
    ],
)
def test_h2o_trimer_multi(levels, component_file, ref_file):

    component_results = load_component_data(component_file)
    ref_data = load_ref_data(ref_file)

    mc = ManyBodyCalculator(mol_h2o_3, [BsseEnum.cp, BsseEnum.nocp, BsseEnum.vmfc], levels, True, False, None)
    nbody_results = mc.analyze(component_results)
    compare_results(nbody_results, ref_data, levels)
