import json
import math
import os

import pytest

from qcmanybody.models import BsseEnum
from qcmanybody.qcengine import run_qcengine
from .util import compare
from .common import mol_h2o_3, specifications

this_dir = os.path.dirname(os.path.abspath(__file__))


@pytest.mark.parametrize(
    "levels, ref_file",
    [
        ({1: "e_mp2", 2: "e_b3lyp", "supersystem": "e_scf"}, "h2o_trimer_multi_ss.json"),
        ({1: "e_mp2", "supersystem": "e_scf"}, "h2o_trimer_multi_ss_2.json"),
    ],
)
def test_h2o_trimer_multi(levels, ref_file):

    nbody_results = run_qcengine(mol_h2o_3, levels, specifications, [BsseEnum.cp, BsseEnum.nocp, BsseEnum.vmfc], True)

    ref_file = os.path.join(this_dir, ref_file)
    with open(ref_file, "r") as f:
        ref_data = json.load(f)

    import pprint
    pprint.pprint(nbody_results)

    compare(nbody_results["ret_energy"], ref_data["CURRENT ENERGY"])

    compare(nbody_results["energy_body_dict"]["1cp"], ref_data["1CP"])
    compare(nbody_results["energy_body_dict"]["1nocp"], ref_data["1NOCP"])
    compare(nbody_results["energy_body_dict"]["1vmfc"], ref_data["1VMFC"])

    if 2 in levels:
        compare(nbody_results["energy_body_dict"]["2cp"], ref_data["2CP"])
        compare(nbody_results["energy_body_dict"]["2nocp"], ref_data["2NOCP"])
        compare(nbody_results["energy_body_dict"]["2vmfc"], ref_data["2VMFC"])

    compare(nbody_results["energy_body_dict"]["3cp"], ref_data["3CP"])
    compare(nbody_results["energy_body_dict"]["3nocp"], ref_data["3NOCP"])
    compare(nbody_results["energy_body_dict"]["3vmfc"], ref_data["3VMFC"])