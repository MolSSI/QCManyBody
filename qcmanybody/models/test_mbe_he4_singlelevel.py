import pprint

import pytest

from qcelemental import constants
from qcelemental.models import Molecule
# v2: from qcelemental.models.procedures_manybody import AtomicSpecification, ManyBodyKeywords, ManyBodyInput
from qcelemental.testing import compare_values

from qcmanybody.models.manybody_v1 import AtomicSpecification, ManyBodyKeywords, ManyBodyInput
from qcmanybody.models.qcng_computer import ManyBodyComputerQCNG, qcvars_to_manybodyproperties

import qcengine as qcng
from addons import using

def skprop(qcvar):
    # qcng: return qcng.procedures.manybody.qcvars_to_manybodyproperties[qcvar]
    return qcvars_to_manybodyproperties[qcvar]


he4_refs_conv = {
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":           -11.530668717083888,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":           -11.522467757090013,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":           -11.522702864080149,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":           -11.522639870651439,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":       0.0,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":       0.008200959993875045,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":       0.007965853003739198,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":       0.008028846432448944,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":           0.008200959993875045,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":          -0.00023510699013584713,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":           6.299342870974556e-05,

        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.530668717083888,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.522851206178828,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.523095269671348,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.523038093664368,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":     0.0,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.007817510905059777,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     0.0075734474125397355,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.007630623419519367,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.007817510905059777,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        -0.00024406349252004134,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         5.717600697963121e-05,

        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.530668717083888,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.52244892169719,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.52268452228489,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.522621528856181,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":     0.0,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.00821979538669737,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     0.007984194798996924,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.00804718822770667,
        "VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.00821979538669737,
        "VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        -0.00023560058770044634,
        "VMFC-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         6.299342870974556e-05,
}

he4_refs_df = {
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":          -11.530751941948,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":          -11.522403579651,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":          -11.522640167467,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":          -11.522576639404,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":      0.0,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":      0.008348362297,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":      0.008111774481,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":      0.008175302544,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":          0.008348362297,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":         -0.000236587816,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":          0.000063528063,

        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":        -11.530751941948,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":        -11.522760073327,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":        -11.523005411447,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":        -11.522948420000,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":    0.0,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":    0.007991868621,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":    0.007746530501,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":    0.007803521948,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":        0.007991868621,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":       -0.000245338120,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":        0.000056991448,

        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY":        -11.530751941948,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY":        -11.522390319401,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY":        -11.522627256726,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY":        -11.522563728663,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":    0.0,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":    0.008361622547,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":    0.008124685222,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":    0.008188213285,
        "VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":        0.008361622547,
        "VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":       -0.000236937325,
        "VMFC-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":        0.000063528063,
    }

sumdict = {
    "4b_all": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cpvmfc": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocpcp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocp_rtd_sio": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocp_sio": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cp_rtd_sio": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cp_sio": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cp_rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_vmfc_rtd": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_vmfc": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",  # TODO remove?
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "3b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_cp_rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_vmfc_rtd": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_vmfc": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY",  # TODO remove?
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "2b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b_nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b_cp_rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b_cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b_vmfc_rtd": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b_vmfc": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY",  # TODO remove?
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "1b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b_nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b_cp_rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b_cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b_vmfc_rtd": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b_vmfc": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY",  # TODO remove?
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
}


@pytest.fixture
def he_tetramer():
    a2 = 2 / constants.bohr2angstroms
    return Molecule(symbols=["He", "He", "He", "He"], fragments=[[0], [1], [2], [3]], geometry=[0, 0, 0, 0, 0, a2, 0, a2, 0, 0, a2, a2])


@pytest.mark.parametrize("program,basis,keywords", [
    pytest.param("cfour", "aug-pvdz", {"frozen_core": False}, id="cfour_conv", marks=using("cfour")),
    pytest.param("gamess", "accd", {"contrl__ispher": 1, "mp2__nacore": 0}, id="gamess_conv", marks=using("gamess")),
    pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "scf__thresh": 1.0e-8, "mp2__freeze": False}, id="nwchem_conv", marks=using("nwchem")),
    pytest.param("psi4", "aug-cc-pvdz", {"e_convergence": 1.e-10, "d_convergence": 1.e-10, "scf_type": "pk", "mp2_type": "conv"}, id="psi4_conv", marks=using("psi4")),
    pytest.param("psi4", "aug-cc-pvdz", {"e_convergence": 1.e-10, "d_convergence": 1.e-10}, id="psi4_df", marks=using("psi4")),
])
@pytest.mark.parametrize("mbe_keywords,anskey,bodykeys,calcinfo_nmbe", [
    pytest.param(
        {"bsse_type": ["nocp", "cp", "vmfc"]},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv],
        65,
        id="4b_all"),
    pytest.param(
        {"bsse_type": ["cp", "vmfc"]},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if ((k.startswith("CP-") and ("TOTAL" not in k)) or (k.startswith("VMFC-")))],
        # TODO: when vmfc active, nocp always available up to max_nbody. cp available if max_nbody=nfr. activate?
        65,
        id="4b_cpvmfc"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "supersystem_ie_only": True},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
        5,
        id="4b_nocp_rtd_sio"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "supersystem_ie_only": True},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
        5,
        id="4b_nocp_sio"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "supersystem_ie_only": True},
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
        9,
        id="4b_cp_rtd_sio"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "supersystem_ie_only": True},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k) and "TOTAL ENERGY" not in k)],
        5,
        id="4b_cp_sio"),
## TODO add vmfc. 3b nmbe=50
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-"))],
        15,
        id="4b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-"))],
        15,
        id="4b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True},
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-"))],
        19,
        id="4b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "TOTAL ENERGY" not in k)],
        15,
        id="4b_cp"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": True},
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-"))],
        65,
        id="4b_vmfc_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False},
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-"))],
        65,  # TODO 61 in reach
        id="4b_vmfc"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 3},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k))],
        14,
        id="3b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 3},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and "4-BODY" not in k)],
        14,
        id="3b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 3},
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "4-BODY" not in k)],
        18,  # bugfix: was 28
        id="3b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 3},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "4-BODY" not in k and "TOTAL ENERGY" not in k)],
        14,
        id="3b_cp"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": True, "max_nbody": 3},
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-") and "4-BODY" not in k)],
        50,
        id="3b_vmfc_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False, "max_nbody": 3},
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-") and "4-BODY" not in k)], # and "TOTAL ENERGY" not in k)],
        50,  # TODO 46 in reach,
        id="3b_vmfc"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 2},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        10,
        id="2b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 2},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        10,
        id="2b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 2},
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        14,
        id="2b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 2},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("4-BODY" not in k) and ("3-BODY" not in k) and "TOTAL ENERGY" not in k)],
        10,
        id="2b_cp"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": True, "max_nbody": 2},
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        22,
        id="2b_vmfc_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False, "max_nbody": 2},
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-") and ("4-BODY" not in k) and ("3-BODY" not in k))], # and "TOTAL ENERGY" not in k)],
        22,  # TODO 18 in reach
        id="2b_vmfc"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 1},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("1-BODY" in k))],
        4,
        id="1b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 1},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("1-BODY" in k))],
        4,  # maybe TODO this could be 0 but rtd hasn't be used to winnow nocp
        id="1b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 1},
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("1-BODY" in k))],
        4,
        id="1b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 1},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("1-BODY" in k) and "TOTAL ENERGY" not in k)],
        0,
        id="1b_cp"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": True, "max_nbody": 1},
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-") and ("1-BODY" in k))],
        4,
        id="1b_vmfc_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False, "max_nbody": 1},
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-") and ("1-BODY" in k))],
        4,  # maybe TODO this could be 0 but rtd hasn't be used to winnow vmfc
        id="1b_vmfc"),
])
def test_nbody_he4_single(program, basis, keywords, mbe_keywords, anskey, bodykeys, calcinfo_nmbe, he_tetramer, request):
    #! MP2/aug-cc-pvDZ many body energies of an arbitrary Helium complex,
    #   addressing 4-body formulas
    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, ...)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    if program == "gamess":
        with pytest.raises(ValueError) as exe:
            # qcng: qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
            ManyBodyComputerQCNG.from_qcschema_ben(mbe_model)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    # qcng: ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    ret = ManyBodyComputerQCNG.from_qcschema_ben(mbe_model)
    print(f"SSSSSSS {request.node.name}")
    # v2: pprint.pprint(ret.model_dump(), width=200)
    pprint.pprint(ret.dict(), width=200)

    _inner = request.node.name.split("[")[1].split("]")[0]
    kwdsln, progln = _inner.split("-")
    refs = he4_refs_df if progln == "psi4_df" else he4_refs_conv
    ans = refs[anskey]
    ref_nmbe = calcinfo_nmbe
    atol = 1.0e-8

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if qcv in bodykeys:
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=atol, label=f"[a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[z] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv in sumdict["4b_all"]:
        skp = skprop(qcv)
        if qcv in sumdict[kwdsln]:
            ref = refs[sumdict[kwdsln][qcv]]
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=atol, label=f"[c] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[d] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[y] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=atol, label=f"[e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=atol, label=f"[g] ret")
    assert ret.properties.calcinfo_nmbe == ref_nmbe, f"{ret.properties.calcinfo_nmbe=} != {ref_nmbe}"
