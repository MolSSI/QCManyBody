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
        ({1: "e_mp2", 2: "e_b3lyp", 3: "e_scf"}, "h2o_trimer_multi.json"),
    ],
)
def test_h2o_trimer_multi(levels, ref_file):

    nbody_results = run_qcengine(mol_h2o_3, levels, specifications, [BsseEnum.cp, BsseEnum.nocp, BsseEnum.vmfc], True)

    ref_file = os.path.join(this_dir, ref_file)
    with open(ref_file, "r") as f:
        ref_data = json.load(f)

    compare(nbody_results["ret_energy"], ref_data["CURRENT ENERGY"])

    compare(nbody_results["energy_body_dict"]["1cp"], ref_data["1CP"])
    compare(nbody_results["energy_body_dict"]["1nocp"], ref_data["1NOCP"])
    compare(nbody_results["energy_body_dict"]["1vmfc"], ref_data["1VMFC"])

    compare(nbody_results["energy_body_dict"]["2cp"], ref_data["2CP"])
    compare(nbody_results["energy_body_dict"]["2nocp"], ref_data["2NOCP"])
    compare(nbody_results["energy_body_dict"]["2vmfc"], ref_data["2VMFC"])

    compare(nbody_results["energy_body_dict"]["3cp"], ref_data["3CP"])
    compare(nbody_results["energy_body_dict"]["3nocp"], ref_data["3NOCP"])
    compare(nbody_results["energy_body_dict"]["3vmfc"], ref_data["3VMFC"])

    ###########################################################

    res = nbody_results["results"]

    compare(res["CP-CORRECTED TOTAL ENERGY"], ref_data["CP-CORRECTED TOTAL ENERGY"])
    compare(res["CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY"], ref_data["CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY"])
    compare(res["CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY"], ref_data["CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY"])
    compare(res["CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY"], ref_data["CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY"])

    compare(res["CP-CORRECTED INTERACTION ENERGY"], ref_data["CP-CORRECTED INTERACTION ENERGY"])
    compare(
        res["CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY"],
        ref_data["CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY"],
    )
    compare(
        res["CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"],
        ref_data["CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"],
    )

    compare(res["CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY"], ref_data["CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY"])
    compare(res["CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY"], ref_data["CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY"])

    ###########################################################

    compare(res["NOCP-CORRECTED TOTAL ENERGY"], ref_data["NOCP-CORRECTED TOTAL ENERGY"])
    compare(res["NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY"], ref_data["NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY"])
    compare(res["NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY"], ref_data["NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY"])
    compare(res["NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY"], ref_data["NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY"])

    compare(res["NOCP-CORRECTED INTERACTION ENERGY"], ref_data["NOCP-CORRECTED INTERACTION ENERGY"])
    compare(
        res["NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY"],
        ref_data["NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY"],
    )
    compare(
        res["NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"],
        ref_data["NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"],
    )

    compare(
        res["NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY"], ref_data["NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY"]
    )
    compare(
        res["NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY"], ref_data["NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY"]
    )

    ###########################################################

    compare(res["VMFC-CORRECTED TOTAL ENERGY"], ref_data["VMFC-CORRECTED TOTAL ENERGY"])
    compare(res["VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY"], ref_data["VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY"])
    compare(res["VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY"], ref_data["VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY"])
    compare(res["VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY"], ref_data["VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY"])

    compare(res["VMFC-CORRECTED INTERACTION ENERGY"], ref_data["VMFC-CORRECTED INTERACTION ENERGY"])
    compare(
        res["VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY"],
        ref_data["VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY"],
    )
    compare(
        res["VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"],
        ref_data["VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"],
    )

    compare(
        res["VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY"], ref_data["VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY"]
    )
    compare(
        res["VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY"], ref_data["VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY"]
    )
