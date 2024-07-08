import pprint

import pytest
from qcelemental import constants
from qcelemental.models import Molecule
from qcelemental.testing import compare, compare_values

from .addons import using, uusing

# these tests use QCEngine compute_procedure and are copied from qcengine/tests/test_qcmanybody.py


@uusing("psi4")
@uusing("qcengine")
@pytest.mark.parametrize(
    "optimizer,bsse_type,sio",
    [
        pytest.param("optking", "none", None, marks=using("optking")),
        pytest.param("genoptking", "none", None, marks=using("optking_genopt")),
        pytest.param("genoptking", "nocp", True, marks=using("optking_genopt")),
        pytest.param("genoptking", "cp", False, marks=using("optking_genopt")),
        pytest.param("geometric", "none", None, marks=using("geometric")),
        pytest.param("gengeometric", "none", None, marks=using("geometric_genopt")),
        pytest.param("gengeometric", "nocp", False, marks=using("geometric_genopt")),
        pytest.param("gengeometric", "cp", True, marks=using("geometric_genopt")),
    ],
)
def test_bsse_opt_hf_trimer(optimizer, bsse_type, sio):

    initial_molecule = Molecule.from_data(
        """
F         -0.04288        2.78905        0.00000
H          0.59079        2.03435        0.00000
--
F         -1.94320       -0.70822        0.00000
H         -1.60642        0.21789       -0.00000
--
F          2.03569       -0.60531       -0.00000
H          1.06527       -0.77673        0.00000
units ang
"""
    )

    at_spec = {
        # schema_name needed for differentiation in genopt
        "schema_name": "qcschema_input",
        "model": {
            "method": "hf",
            "basis": "6-31g",
        },
        "keywords": {
            "scf_type": "df",
        },
    }

    mbe_spec = {
        # schema_name needed for differentiation in genopt
        "schema_name": "qcschema_manybodyspecification",
        "specification": {
            "model": {
                "method": "hf",
                "basis": "6-31g",
            },
            "driver": "energy",
            "program": "psi4",
            "keywords": {},
            "protocols": {
                "stdout": False,
            },
            "extras": {
                "psiapi": True,
            },
        },
        "keywords": {
            "bsse_type": bsse_type,
            "supersystem_ie_only": sio,
        },
        "driver": "energy",
        "protocols": {
            "component_results": "all",
        },
    }

    opt_data = {
        "initial_molecule": initial_molecule,
        "input_specification": at_spec if (bsse_type == "none") else mbe_spec,
        "keywords": {
            "program": "psi4",
            "g_convergence": "nwchem_loose",
        },
        "protocols": {
            "trajectory": "initial_and_final",
        },
    }
    # from qcmanybody.models.generalized_optimization import GeneralizedOptimizationInput
    # opt_data = GeneralizedOptimizationInput(**opt_data)

    import qcengine as qcng
    ret = qcng.compute_procedure(opt_data, optimizer, raise_error=True)

    print("FFFFFFFFFF")
    pprint.pprint(ret.dict(), width=200)

    r_fh_hb_xptd = {
        "none": 2.18 / constants.bohr2angstroms,
        "nocp": 2.18 / constants.bohr2angstroms,
        "cp": 2.27 / constants.bohr2angstroms,
    }[bsse_type]
    r_fh_computed = ret.final_molecule.measure([1, 3])
    assert (
        pytest.approx(r_fh_computed, 1.0e-2) == r_fh_hb_xptd
    ), f"hydrogen bond length computed ({r_fh_computed}) != expected ({r_fh_hb_xptd})"
    assert (
        len(ret.trajectory) == 2
    ), f"trajectory protocol did not take. len(ret.trajectory)={len(ret.trajectory)} != 2 (initial_and_final)"
    if bsse_type != "none":
        xptd_nmbe = {
            ("nocp", False): 7,
            ("nocp", True): 4,
            ("cp", False): 10,
            ("cp", True): 7,
        }[(bsse_type, sio)]
        assert xptd_nmbe == len(ret.trajectory[-1].component_results), f"mbe protocol did not take"
        assert (
            ret.trajectory[-1].component_results['["(auto)", [1, 2, 3], [1, 2, 3]]'].stdout is None
        ), f"atomic protocol did not take"


@uusing("qcengine")
@pytest.mark.parametrize("bsse_type", ["mbe", "ssfc"])  # aka nocp, cp
@pytest.mark.parametrize(
    "qcprog,qc_keywords",
    [
        pytest.param("psi4", {}, marks=using("psi4")),
        pytest.param("cfour", {}, marks=using("cfour")),
        pytest.param("nwchem", {}, marks=using("nwchem")),
    ],
)
@pytest.mark.parametrize(
    "optimizer,opt_keywords",
    [
        pytest.param("optking", {}, marks=using("optking_genopt")),
        pytest.param(
            "geometric", {"convergence_set": "interfrag_tight", "maxiter": 15}, marks=using("geometric_genopt")
        ),
    ],
)
def test_bsse_opt_lif_dimer(optimizer, opt_keywords, bsse_type, qcprog, qc_keywords):
    # in geomeTRIC: test_lif_bsse

    lif = {
        "symbols": ["Li", "F"],
        "geometry": [0, 0, 0, 0, 0, 3],
        "fragments": [[0], [1]],
        "fragment_charges": [+1, -1],
    }

    mbe_spec = {
        "schema_name": "qcschema_manybodyspecification",
        "specification": {
            "model": {
                "method": "hf",
                "basis": "6-31G",
            },
            "driver": "energy",
            "program": qcprog,
            "keywords": qc_keywords,
            "protocols": {
                "stdout": False,
            },
        },
        "keywords": {
            "bsse_type": bsse_type,
            "supersystem_ie_only": True,
        },
        "protocols": {
            "component_results": "all",
        },
        "driver": "energy",
    }

    opt_data = {
        "initial_molecule": lif,
        "input_specification": mbe_spec,
        "keywords": {
            "program": "nonsense",
            **opt_keywords,
        },
        "protocols": {
            "trajectory": "final",
        },
    }

    import qcengine as qcng
    ret = qcng.compute_procedure(opt_data, "gen" + optimizer, raise_error=True)

    # printing will show up if job fails
    pprint.pprint(ret.dict(), width=200)

    assert ret.success

    assert len(ret.trajectory[0].component_properties) == (5 if bsse_type == "ssfc" else 3)

    assert ret.provenance.creator.lower() == optimizer
    assert ret.trajectory[0].provenance.creator == "QCManyBody"
    atres = list(ret.trajectory[0].component_results.values())[0]
    assert atres.provenance.creator.lower() == qcprog

    Rlif = ret.final_molecule.measure([0, 1])
    Rref = 3.016 if bsse_type == "ssfc" else 2.969
    assert compare_values(Rlif, Rref, "bond length", atol=1.0e-3)
