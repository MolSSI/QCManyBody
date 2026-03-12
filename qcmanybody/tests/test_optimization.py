import pprint

import pytest
from qcelemental import constants
from qcelemental.testing import compare_values

from .addons import using, uusing

# these tests use QCEngine compute_procedure and are copied from qcengine/tests/test_qcmanybody.py


@uusing("psi4")
@uusing("qcengine")
@pytest.mark.parametrize("schver", [1, 2])
@pytest.mark.parametrize(
    "optimizer,bsse_type,sio",
    [
        pytest.param("optking", "nil", None, marks=using("optking")),
        pytest.param("optking", "nocp", True, marks=using("optking_v2")),
        pytest.param("optking", "cp", False, marks=using("optking_v2")),
        pytest.param("geometric", "nil", None, marks=using("geometric")),
        pytest.param("gengeometric", "nil", None, marks=using("geometric_genopt")),
        pytest.param("gengeometric", "nocp", False, marks=using("geometric_genopt")),
        pytest.param("gengeometric", "cp", True, marks=using("geometric_genopt")),
    ],
)
def test_bsse_opt_hf_trimer(optimizer, bsse_type, sio, schver):
    if bsse_type != "nil" and schver == 1:
        pytest.skip("ManyBody Optimization is only available for QCSchema v2. The experimental v1 GeneralizedOptimization is retired.")

    if schver == 1:
        from qcelemental.models.v1 import Molecule
        subptcl = "component_results"
        subres = "component_results"
    elif schver == 2:
        from qcelemental.models.v2 import Molecule
        subptcl = "cluster_results"
        subres = "cluster_results"

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
        "schema_name": "qcschema_atomic_specification" if schver == 2 else "qcschema_input",
        "program": "psi4",
        "driver": "gradient",
        "model": {
            "method": "hf",
            "basis": "6-31g",
        },
        "keywords": {
            "scf_type": "df",
        },
    }
    if schver == 1:
        at_spec.pop("program")
        at_spec.pop("driver")

    mbe_spec = {
        # schema_name needed for differentiation in genopt
        "schema_name": "qcschema_many_body_specification" if schver == 2 else "qcschema_manybodyspecification",
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
            subptcl: "all",
        },
    }

    if schver == 2:
        opt_data = {
            "initial_molecule": initial_molecule,
            "specification": {
                "specification": at_spec if (bsse_type == "nil") else mbe_spec,
                "keywords": {
                    "g_convergence": "nwchem_loose",
                },
                "protocols": {
                    "trajectory_results": "initial_and_final",
                },
            },
        }
    else:
        opt_data = {
            "initial_molecule": initial_molecule,
            "input_specification": at_spec if (bsse_type == "nil") else mbe_spec,
            "keywords": {
                "program": "psi4",
                "g_convergence": "nwchem_loose",
            },
            "protocols": {
                "trajectory": "initial_and_final",
            },
        }

    import qcengine as qcng
    ret = qcng.compute_procedure(opt_data, optimizer, raise_error=True, return_version=schver)

    print("FFFFFFFFFF")
    pprint.pprint(ret.model_dump(), width=200)

    r_fh_hb_xptd = {
        "nil": 2.18 / constants.bohr2angstroms,
        "nocp": 2.18 / constants.bohr2angstroms,
        "cp": 2.27 / constants.bohr2angstroms,
    }[bsse_type]
    r_fh_computed = ret.final_molecule.measure([1, 3])
    assert (
        pytest.approx(r_fh_computed, 1.0e-2) == r_fh_hb_xptd
    ), f"hydrogen bond length computed ({r_fh_computed}) != expected ({r_fh_hb_xptd})"
    ret_traj = ret.trajectory_results if schver == 2 else ret.trajectory
    assert (
        len(ret_traj) == 2
    ), f"trajectory protocol did not take. len(ret.trajectory)={len(ret_traj)} != 2 (initial_and_final)"
    if bsse_type != "nil":
        xptd_nmbe = {
            ("nocp", False): 7,
            ("nocp", True): 4,
            ("cp", False): 10,
            ("cp", True): 7,
        }[(bsse_type, sio)]
        ret_last_subres = getattr(ret_traj[-1], subres)
        assert xptd_nmbe == len(ret_last_subres), f"mbe protocol did not take"
        assert ret_last_subres['["(auto)", [1, 2, 3], [1, 2, 3]]'].stdout is None, f"atomic protocol did not take"


@uusing("qcengine")
@pytest.mark.parametrize("schver", [2])  # GeneralizedOptimization retired (this was the v1 precursor)
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
        pytest.param("optking", {}, marks=using("optking_v2")),
        pytest.param(
            "geometric", {"convergence_set": "interfrag_tight", "maxiter": 15}, marks=using("geometric_genopt")
        ),
    ],
)
def test_bsse_opt_lif_dimer(optimizer, opt_keywords, bsse_type, qcprog, qc_keywords, schver):
    # in geomeTRIC: test_lif_bsse

    lif = {
        "symbols": ["Li", "F"],
        "geometry": [0, 0, 0, 0, 0, 3],
        "fragments": [[0], [1]],
        "fragment_charges": [+1, -1],
    }

    mbe_subprop = "component_properties" if schver == 1 else "cluster_properties"
    mbe_subptcl = "component_results" if schver == 1 else "cluster_results"
    opt_subptcl = "trajectory" if schver == 1 else "trajectory_results"

    mbe_spec = {
        "schema_name": "qcschema_many_body_specification" if schver == 2 else "qcschema_manybodyspecification",
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
            mbe_subptcl: "all",
        },
        "driver": "energy",
    }

    if schver == 2:
      opt_data = {
        "initial_molecule": lif,
        "specification": {
            "specification": at_spec if (bsse_type == "nil") else mbe_spec,
            "keywords": {
                **opt_keywords,
            },
            "protocols": {
                opt_subptcl: "final",
            },
        },
      }
    else:
      opt_data = {
        "initial_molecule": lif,
        "input_specification": at_spec if (bsse_type == "nil") else mbe_spec,
        "keywords": {
            "program": "nonsense",
            **opt_keywords,
        },
        "protocols": {
            opt_subptcl: "final",
        },
      }

    import qcengine as qcng
    ret = qcng.compute_procedure(opt_data, optimizer, raise_error=True, return_version=schver)

    # printing will show up if job fails
    pprint.pprint(ret.model_dump(), width=200)

    assert ret.success

    ret_traj = ret.trajectory_results if schver == 2 else ret.trajectory
    assert len(getattr(ret_traj[0], mbe_subprop)) == (5 if bsse_type == "ssfc" else 3)

    assert ret.provenance.creator.lower() == optimizer
    assert ret_traj[0].provenance.creator == "QCManyBody"
    atres = list(getattr(ret_traj[0], mbe_subptcl).values())[0]
    assert atres.provenance.creator.lower() == qcprog

    Rlif = ret.final_molecule.measure([0, 1])
    Rref = 3.016 if bsse_type == "ssfc" else 2.969
    assert compare_values(Rlif, Rref, "bond length", atol=1.0e-3)
