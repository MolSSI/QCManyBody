import re
import copy
import pprint

import pytest
import numpy as np

from qcelemental.models import Molecule
from qcelemental.testing import compare_values, compare_recursive

from qcmanybody.models.manybody_pydv1 import AtomicSpecification, ManyBodyKeywords, ManyBodyInput
from qcmanybody.models.qcng_computer import ManyBodyComputerQCNG, qcvars_to_manybodyproperties

#import qcengine as qcng
from .addons import using

def skprop(qcvar):
    # qcng: return qcng.procedures.manybody.qcvars_to_manybodyproperties[qcvar]
    return qcvars_to_manybodyproperties[qcvar]


@pytest.fixture(scope="function")
def mbe_data_grad_dtz():
    # note that spherical/cartesian irrelevant for He & 6-31G, and fc/ae irrelevant for He
#    c4_kwds = {}
#    gms_kwds = {"basis__ngauss": 6, "ccinp__ncore": 0, "ccinp__iconv": 9, "scf__conv": 9}
#    nwc_kwds = {"scf__thresh": 1.0e-8, "ccsd__thresh": 1.e-8}
    p4_kwds = {"scf_type": "pk", "mp2_type": "conv"}

    protocols = {"stdout": False}
    return {
        "specification": {
            "specification": {
                "p4-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "cc-pvdz",
                    },
                    "driver": "energy",
                    "program": "psi4",
                    "keywords": p4_kwds,
                    "protocols": protocols,
                },
                "p4-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "cc-pvdz",
                    },
                    "driver": "energy",
                    "program": "psi4",
                    "keywords": p4_kwds,
                    "protocols": protocols,
                },
                "p4-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "cc-pvdz",
                    },
                    "driver": "energy",
                    "program": "psi4",
                    "keywords": p4_kwds,
                    "protocols": protocols,
                },
            },
            "keywords": None,
            "driver": "energy",
        },
        "molecule": None,
    }


het4_refs_conv_grad_dtz = {
    # 4: hf/cc-pvdz ; copied from psi4 - not computed from pieces
    "4": {
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":       -203.05040290887501,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":       -203.35196718523923,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":       -203.32671800908514,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":       -203.32531316984924,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   -0.3015642763642177,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   -0.27631510021012673,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   -0.2749102609742238,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       -0.3015642763642177,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        0.025249176154090947,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":        0.001404839235902955,

        "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY": np.array([
             [-7.40362091e-02,  2.84130342e-03,  1.00812419e-02],
             [ 7.57644607e-02, -1.48383966e-04, -1.73598306e-03],
             [ 2.60812740e-02,  4.40153410e-04,  2.02181413e-03],
             [-8.52081479e-03, -2.26666899e-02, -1.10830218e-01],
             [-6.03761361e-03, -1.54080944e-02,  1.00510650e-01],
             [-1.32510973e-02,  3.49417115e-02, -4.75050095e-05]]),
    }
}


#EE through 3b
#        n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#                   [Eh]                    [Eh]                  [kcal/mol]            [Eh]                  [kcal/mol]
#             1     -203.434092764220                   nan                   nan        0.000000000000        0.000000000000
#             2     -203.338005606277                   nan                   nan        0.096087157943       60.295601917486
#       RTN   3     -203.326023705227                   nan                   nan        0.011981901050        7.518756422798
#      FULL   4        N/A                   N/A                   N/A                   N/A                   N/A             


# only here for keys
he4_refs_conv = {
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":           None,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":           None,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":           None,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":           None,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":     None,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     None,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     None,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     None,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         None,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":         None,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         None,

        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         None,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         None,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         None,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         None,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":   None,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   None,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   None,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   None,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       None,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":       None,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":       None,

        "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY":       None,

        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         None,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         None,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         None,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         None,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":   None,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   None,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   None,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   None,
        "VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       None,
        "VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":       None,
        "VMFC-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":       None,
}


sumstr = {
    "4b_nocp_rtd": {
        "4": r"""
   ==> N-Body: Non-Counterpoise Corrected \(NoCP\) energies <==

^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+1\s+-203.0504029\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+2\s+-203.3519671\d+       -0.3015642\d+     -189.23\d+       -0.3015642\d+     -189.23\d+
^\s+3\s+-203.3267180\d+       -0.2763151\d+     -173.39\d+        0.0252491\d+       15.84\d+
^\s+FULL/RTN\s+4     -203.3253131\d+       -0.2749102\d+     -172.50\d+        0.0014048\d+        0.88\d+
""",
#   "cp4b_tot": {
#       "121": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+2\s+-11.4710585\d+        0.0095899\d+        6.01\d+        0.0095899\d+        6.01\d+
#^\s+3\s+-11.4713246\d+        0.0093239\d+        5.85\d+       -0.0002660\d+       -0.16\d+
#^\s+FULL/RTN\s+4\s+-11.4712722\d+        0.0093763\d+        5.88\d+        0.0000523\d+        0.03\d+
#""",
#        "22": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+2\s+-11.4707059\d+        0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
#^\s+3\s+-11.4709719\d+        0.0096765\d+        6.07\d+       -0.0002660\d+       -0.16\d+
#^\s+FULL/RTN\s+4\s+-11.4709134\d+        0.0097351\d+        6.10\d+        0.0000585\d+        0.03\d+
#""",
#        "ss22": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+2\s+-11.4707059\d+        0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
#^\s+SS\s+3\s+N/A\s+N/A\s+N/A\s+N/A\s+N/A\s*
#^\s+SS/FULL/RTN\s+4      -11.4707760\d+        0.0098725\d+        6.19\d+       -0.0000700\d+       -0.04\d+
#""",
#   },
#   "cp3b_tot": {
#       "121": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+2\s+-11.4710585\d+        0.0095899\d+        6.01\d+        0.0095899\d+        6.01\d+
#^\s+RTN\s+3\s+-11.4713246\d+        0.0093239\d+        5.85\d+       -0.0002660\d+       -0.16\d+
#^\s+FULL\s+4\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
#""",
#        "22": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+2\s+-11.4707059\d+        0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
#^\s+RTN\s+3\s+-11.4709719\d+        0.0096765\d+        6.07\d+       -0.0002660\d+       -0.16\d+
#^\s+FULL\s+4\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
#""",
#   },
#   "cp3b_ie": {
#       "121": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+N/A\s+0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+2\s+N/A\s+0.0095899\d+        6.01\d+        0.0095899\d+        6.01\d+
#^\s+RTN\s+3\s+N/A\s+0.0093239\d+        5.85\d+       -0.0002660\d+       -0.16\d+
#^\s+FULL\s+4        N/A                   N/A                   N/A                   N/A                   N/A\s*
#""",
#        "22": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+N/A\s+0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+2\s+N/A\s+0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
#^\s+RTN\s+3\s+N/A\s+0.0096765\d+        6.07\d+       -0.0002660\d+       -0.16\d+
#^\s+FULL\s+4        N/A                   N/A                   N/A                   N/A                   N/A\s*
#""",
#   },
#   "cp2b_tot": {
#       "121": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+RTN\s+2\s+-11.4710585\d+        0.0095899\d+        6.01\d+        0.0095899\d+        6.01\d+
#^\s+3\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
#^\s+FULL\s+4\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
#""",
#        "22": r"""
#   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
#
#^\s+n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
#^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
#^\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
#^\s+RTN\s+2\s+-11.4707059\d+        0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
#^\s+3\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
#^\s+FULL\s+4\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
#""",
   },
}

sumdict_grad = {
    "4b_all": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",

        "NOCP-CORRECTED TOTAL GRADIENT": "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY",
        #"NOCP-CORRECTED INTERACTION GRADIENT": "NOCP-CORRECTED INTERACTION GRADIENT THROUGH 4-BODY",
        "CP-CORRECTED TOTAL GRADIENT": "CP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY",
        #"CP-CORRECTED INTERACTION GRADIENT": "CP-CORRECTED INTERACTION GRADIENT THROUGH 4-BODY",
        "VMFC-CORRECTED TOTAL GRADIENT": "VMFC-CORRECTED TOTAL GRADIENT THROUGH 4-BODY",
        #"VMFC-CORRECTED INTERACTION GRADIENT": "VMFC-CORRECTED INTERACTION GRADIENT THROUGH 4-BODY",
    },
    "4b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",

        "NOCP-CORRECTED TOTAL GRADIENT": "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY",
        #"NOCP-CORRECTED INTERACTION GRADIENT": "NOCP-CORRECTED INTERACTION GRADIENT THROUGH 4-BODY",
    },
}


@pytest.fixture
def het_tetramer():
    return Molecule(
        symbols=["F", "H", "F", "H", "H", "He"], 
        fragments=[[0], [1, 2], [3, 4], [5]], 
        fragment_charges=[-1, 0, 0, 0],
        geometry=[-2,  0,  0, 0,  0,  0, 4,  0,  0, 0,  3,  0, 0,  3,  2, 0, -3,  0],
    )


@pytest.mark.parametrize("levels", [
    # pattern 4
    pytest.param({4: "p4-hf"}, id="4-psi4_pure", marks=using("psi4")),

])
@pytest.mark.parametrize("mbe_keywords,anskeyE,anskeyG,bodykeys,outstrs,calcinfo_nmbe", [
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY",
        [k for k in he4_refs_conv if k.startswith("NOCP-")],
        ["4b_nocp_rtd"],
        {"4": 15,},
        #All,4b QCVariable: 1_((1, 2, 3, 4), (1, 2, 3, 4)) -203.32531316984924
        id="4b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "embedding_charges": {1: [-1.0], 2: [0.5, -0.5], 3: [-0.5, 0.5], 4: [0]}},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY",
        [k for k in he4_refs_conv if k.startswith("NOCP-")],
        None, #["4b_nocp_rtd"],
        {"4": 15,},
        #All,4b QCVariable: 1_((1, 2, 3, 4), (1, 2, 3, 4)) -203.32531316984924
        id="4b_nocp_rtd_ee"),
])
def test_nbody_het4_grad(levels, mbe_keywords, anskeyE, anskeyG, bodykeys, outstrs, calcinfo_nmbe, het_tetramer, request, mbe_data_grad_dtz):
    _inner = request.node.name.split("[")[1].split("]")[0]
    kwdsln, pattern, progln = _inner.split("-")

    levels = copy.deepcopy(levels)

    mbe_keywords = ManyBodyKeywords(levels=levels, **mbe_keywords)
    mbe_data_grad_dtz["molecule"] = het_tetramer
    mbe_data_grad_dtz["specification"]["driver"] = "gradient"
    mbe_data_grad_dtz["specification"]["keywords"] = mbe_keywords
    mbe_model = ManyBodyInput(**mbe_data_grad_dtz)

    # qcng: ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    ret = ManyBodyComputerQCNG.from_qcschema_ben(mbe_model)
    print(f"MMMMMMM {request.node.name}")
    pprint.pprint(ret.dict(), width=200)

    refs = het4_refs_conv_grad_dtz[pattern]
    ansE = refs[anskeyE]
    ansG = refs[anskeyG]
    ref_nmbe = calcinfo_nmbe[pattern]
    ref_bodykeys = bodykeys[pattern] if pattern in bodykeys else bodykeys
    ref_sumdict = sumdict_grad[kwdsln][pattern] if pattern in sumdict_grad[kwdsln] else sumdict_grad[kwdsln]
    atol = 2.5e-8

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if qcv in ref_bodykeys:
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=atol, label=f"[a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[z] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv in sumdict_grad["4b_all"]:
        skp = skprop(qcv)
        if qcv in ref_sumdict:
            ref = refs[ref_sumdict[qcv]]
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=atol, label=f"[c] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[d] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[y] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CURRENT ENERGY": ansE,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=atol, label=f"[e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[f] skprop {skp}")

    for qcv, ref in {
        "CURRENT GRADIENT": ansG,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=atol, label=f"[e] G qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[f] G skprop {skp}")
    assert compare_values(ansG, ret.return_result, atol=atol, label=f"[g] G ret")

    assert ret.properties.calcinfo_nmbe == ref_nmbe, f"{ret.properties.calcinfo_nmbe=} != {ref_nmbe}"

    if outstrs:
        for stdoutkey in outstrs:
            assert re.search(sumstr[stdoutkey][pattern], ret.stdout, re.MULTILINE), f"[j] N-Body pattern not found: {sumstr[stdoutkey][pattern]}"
