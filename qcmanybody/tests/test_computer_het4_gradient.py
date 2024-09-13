import os
import pprint
import re

import numpy as np
import pytest
from qcelemental.models import Molecule
from qcelemental.testing import compare_values

from qcmanybody.computer import ManyBodyComputer
from qcmanybody.models import ManyBodyInput, ManyBodyKeywords, ManyBodyResultProperties
from qcmanybody.utils import translate_qcvariables

from .addons import uusing


@pytest.fixture(scope="function")
def mbe_data_grad_dtz():
    # note that spherical/cartesian irrelevant for He & 6-31G, and fc/ae irrelevant for He
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
    "4b_nocp_rtd": {
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
    },
    "4b_nocp_rtd_ee": {
        # through 3b from psi4. 4b from qcmb!
            # TODO IE should be present"
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":       -203.4340927642202,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":       -203.33800560627736,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":       -203.32602370522727,
        #"NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":    0.09608715794283285,
        #"NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":    0.10806905899292474,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":        0.09608715794283285,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        0.011981901050091892,
        #"NOCP-CORRECTED TOTAL GRADIENT THROUGH 3-BODY": np.array([  # fine but "wrongly present"
        #    [-7.36025774e-02,  3.36991350e-03,  1.00686673e-02],
        #    [ 7.55655033e-02, -5.37981728e-04, -2.19873195e-03],
        #    [ 2.63105970e-02,  4.12475495e-04,  2.10686911e-03],
        #    [-8.78129292e-03, -2.22190993e-02, -1.10826802e-01],
        #    [-6.14121135e-03, -1.53249058e-02,  1.00569730e-01],
        #    [-1.32579782e-02,  3.44008392e-02,  3.46216252e-05]]),

        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":       -203.325313169849,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":        0.000710535378,
        "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY": np.array([
            [-7.40362107e-02,  2.84130177e-03,  1.00812417e-02],
            [ 7.57644630e-02, -1.48382170e-04, -1.73598308e-03],
            [ 2.60812737e-02,  4.40153586e-04,  2.02181449e-03],
            [-8.52081509e-03, -2.26666902e-02, -1.10830218e-01],
            [-6.03761335e-03, -1.54080947e-02,  1.00510650e-01],
            [-1.32510976e-02,  3.49417117e-02, -4.75050664e-05]]),
    },
}



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
    "4b_nocp_rtd": r"""
   ==> N-Body: Non-Counterpoise Corrected \(NoCP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-203.0504029\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+§A\s+2\s+-203.3519671\d+       -0.3015642\d+     -189.23\d+       -0.3015642\d+     -189.23\d+
^\s+§A\s+3\s+-203.3267180\d+       -0.2763151\d+     -173.39\d+        0.0252491\d+       15.84\d+
^\s+FULL/RTN\s+§A\s+4     -203.3253131\d+       -0.2749102\d+     -172.50\d+        0.0014048\d+        0.88\d+
""",
    "4b_nocp_rtd_ee": r"""
^\s+==> N-Body: Non-Counterpoise Corrected \(NoCP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-203.4340927\d+\s+N/A\s+N/A\s+0.0000000\d+        0.00\d+
^\s+§A\s+2\s+-203.3380056\d+\s+N/A\s+N/A\s+0.0960871\d+       60.29\d+
^\s+§A\s+3\s+-203.3260237\d+\s+N/A\s+N/A\s+0.0119819\d+        7.51\d+
^\s+FULL/RTN\s+§A\s+4\s+-203.3253131\d+\s+N/A\s+N/A\s+0.0007105\d+\s+0.44\d+
""",
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
    "4b_nocp_rtd_ee": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        #"NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",

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


@uusing("psi4")
@pytest.mark.parametrize("mbe_keywords,anskeyE,anskeyG,bodykeys,outstrs,calcinfo_nmbe", [
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "levels": {4: "p4-hf"}},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY",
        [k for k in he4_refs_conv if k.startswith("NOCP-")],
        ["4b_nocp_rtd"],
        15,
        #All,4b QCVariable: 1_((1, 2, 3, 4), (1, 2, 3, 4)) -203.32531316984924
        id="4b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "levels": {4: "p4-hf"}, "embedding_charges": {1: [-1.0], 2: [0.5, -0.5], 3: [-0.5, 0.5], 4: [0]}},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED TOTAL GRADIENT THROUGH 4-BODY",
        [k for k in he4_refs_conv if k.startswith("NOCP-")],
        ["4b_nocp_rtd_ee"],
        15,
        id="4b_nocp_rtd_ee"),
])
def test_nbody_het4_grad(mbe_keywords, anskeyE, anskeyG, bodykeys, outstrs, calcinfo_nmbe, het_tetramer, request, mbe_data_grad_dtz, monkeypatch):
    _inner = request.node.name.split("[")[1].split("]")[0]
    kwdsln, pattern, progln = _inner, "", "psi4"
    monkeypatch.setenv("QCMANYBODY_EMBEDDING_CHARGES", "1")

    mbe_keywords = ManyBodyKeywords(**mbe_keywords)
    mbe_data_grad_dtz["molecule"] = het_tetramer
    mbe_data_grad_dtz["specification"]["driver"] = "gradient"
    mbe_data_grad_dtz["specification"]["keywords"] = mbe_keywords
    mbe_model = ManyBodyInput(**mbe_data_grad_dtz)

    # qcng: ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    ret = ManyBodyComputer.from_manybodyinput(mbe_model)
    print(f"MMMMMMM {request.node.name}")
    pprint.pprint(ret.dict(), width=200)

    refs = het4_refs_conv_grad_dtz[kwdsln]
    ansE = refs[anskeyE]
    ansG = refs[anskeyG]
    ref_nmbe = calcinfo_nmbe
    ref_bodykeys = bodykeys
    ref_sumdict = sumdict_grad[kwdsln]
    atol = 2.5e-8

    # don't want QCVariables stashed in extras, but prepare the qcvars translation, and check it
    assert ret.extras == {}, f"[w] extras wrongly present: {ret.extras.keys()}"
    qcvars = translate_qcvariables(ret.properties.dict())

    skprop = ManyBodyResultProperties.to_qcvariables(reverse=True)

    for qcv, ref in refs.items():
        skp = skprop[qcv]
        if qcv in ref_bodykeys:
            assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[b] skprop {skp}")
        else:
            assert qcv not in qcvars, f"[z] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv in sumdict_grad["4b_all"]:
        skp = skprop[qcv]
        if qcv in ref_sumdict:
            ref = refs[ref_sumdict[qcv]]
            assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[c] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[d] skprop {skp}")
        else:
            assert qcv not in qcvars, f"[y] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CURRENT ENERGY": ansE,
    }.items():
        skp = skprop[qcv]
        assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[f] skprop {skp}")

    for qcv, ref in {
        "CURRENT GRADIENT": ansG,
    }.items():
        skp = skprop[qcv]
        assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[e] G qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[f] G skprop {skp}")
    assert compare_values(ansG, ret.return_result, atol=atol, label=f"[g] G ret")

    assert ret.properties.calcinfo_nmbe == ref_nmbe, f"{ret.properties.calcinfo_nmbe=} != {ref_nmbe}"

    if outstrs:
        for stdoutkey in outstrs:
            assert re.search(sumstr[stdoutkey], ret.stdout, re.MULTILINE), f"[j] N-Body pattern not found: {sumstr[stdoutkey]}"


@pytest.mark.parametrize("mbe_keywords,errmsg", [
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "levels": {4: "p4-hf"}, "embedding_charges": {1: [-1.0], 2: [0.5, -0.5], 3: [-0.5, 0.5], 4: [0]}},
        "Embedding charges for EE-MBE are still in testing",
        id="4b_nocp_rtd_ee_error"),
])
def test_nbody_ee_error(mbe_keywords, errmsg, het_tetramer, mbe_data_grad_dtz):

    mbe_keywords = ManyBodyKeywords(**mbe_keywords)
    mbe_data_grad_dtz["molecule"] = het_tetramer
    mbe_data_grad_dtz["specification"]["driver"] = "gradient"
    mbe_data_grad_dtz["specification"]["keywords"] = mbe_keywords
    mbe_model = ManyBodyInput(**mbe_data_grad_dtz)

    with pytest.raises(ValueError) as e:
        ManyBodyComputer.from_manybodyinput(mbe_model)

    assert errmsg in str(e.value), e.value


def test_fragmentless_mol(mbe_data_grad_dtz):
    het_tetramer_fragmentless = Molecule(
        symbols=["F", "H", "F", "H", "H", "He"],
        geometry=[-2,  0,  0, 0,  0,  0, 4,  0,  0, 0,  3,  0, 0,  3,  2, 0, -3,  0],
    )

    mbe_data_grad_dtz["molecule"] = het_tetramer_fragmentless
    mbe_data_grad_dtz["specification"]["keywords"] = {}
    mbe_model = ManyBodyInput(**mbe_data_grad_dtz)

    with pytest.raises(ValueError) as e:
        ManyBodyComputer.from_manybodyinput(mbe_model)

    assert "fragmentation has not been specified" in str(e.value), e.value
