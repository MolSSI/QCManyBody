import os

import pytest

from qcmanybody import ManyBodyCore
from qcmanybody.models import BsseEnum

from .common import mol_h2o_3
from .utils import compare_results, load_component_data, load_ref_data

this_dir = os.path.dirname(os.path.abspath(__file__))


@pytest.mark.parametrize(
    "levels, component_file, ref_file",
    [
        # fmt: off
        # Grad reference file includes energies
        ({1: "e_mp2", 2: "e_b3lyp", "supersystem": "e_scf"}, "h2o_trimer_multi_ss_energy_1",   "h2o_trimer_multi_ss_gradient_1"),
        ({1: "e_mp2",               "supersystem": "e_scf"}, "h2o_trimer_multi_ss_energy_2",   "h2o_trimer_multi_ss_gradient_2"),
        ({1: "g_mp2", 2: "g_b3lyp", "supersystem": "g_scf"}, "h2o_trimer_multi_ss_gradient_1", "h2o_trimer_multi_ss_gradient_1"),
        ({1: "g_mp2",               "supersystem": "g_scf"}, "h2o_trimer_multi_ss_gradient_2", "h2o_trimer_multi_ss_gradient_2"),
        # ({1: "h_scf_atz", 2: "h_scf_adz", 3: "h_scf"}, "h2o_trimer_multi_ss_hessian_1", "h2o_trimer_multi_ss_hessian_1"),
        # ({1: "h_scf_atz", 2: "h_scf_adz"            }, "h2o_trimer_multi_ss_hessian_2", "h2o_trimer_multi_ss_hessian_2"),
        # fmt: on
    ],
)
def test_h2o_trimer_multi_ss(levels, component_file, ref_file):
    component_results = load_component_data(component_file)
    ref_data = load_ref_data(ref_file)

    mc = ManyBodyCore(mol_h2o_3, [BsseEnum.cp, BsseEnum.nocp, BsseEnum.vmfc], levels,
        return_total_data=True, supersystem_ie_only=False, embedding_charges=None)
    nbody_results = mc.analyze(component_results)
    compare_results(nbody_results, ref_data, levels)
