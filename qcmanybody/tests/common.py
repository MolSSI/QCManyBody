from qcelemental.models import Molecule

mol_ne2 = Molecule(
    symbols=["Ne", "Ne"], fragments=[[0], [1]], geometry=[0.0, 0.0, -2.834589188186742, 0.0, 0.0, 2.834589188186742]
)

mol_h2o_3 = Molecule.from_data(
    """
O      -2.76373224  -1.24377706  -0.15444566
H      -1.12357791  -2.06227970  -0.05243799
H      -3.80792362  -2.08705525   1.06090407
--
O       2.46924614  -1.75437739  -0.17092884
H       3.76368260  -2.21425403   1.00846104
H       2.30598330   0.07098445  -0.03942473
--
O       0.29127930   3.00875625   0.20308515
H      -1.21253048   1.95820900   0.10303324
H       0.10002049   4.24958115  -1.10222079
units bohr
"""
)

mol_h2o_3_dict = {k: v for k, v in mol_h2o_3.dict().items() if k in ["symbols", "geometry", "fragments"]}

specifications = {
    "e_scf": {
        "program": "psi4",
        "specification": {
            "driver": "energy",
            "model": {"method": "hf", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "e_mp2": {
        "program": "psi4",
        "specification": {
            "driver": "energy",
            "model": {"method": "mp2", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "e_b3lyp": {
        "program": "psi4",
        "specification": {
            "driver": "energy",
            "model": {"method": "b3lyp", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "g_scf": {
        "program": "psi4",
        "specification": {
            "driver": "gradient",
            "model": {"method": "hf", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "g_mp2": {
        "program": "psi4",
        "specification": {
            "driver": "gradient",
            "model": {"method": "mp2", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "g_b3lyp": {
        "program": "psi4",
        "specification": {
            "driver": "gradient",
            "model": {"method": "b3lyp", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "h_scf": {
        "program": "psi4",
        "specification": {
            "driver": "hessian",
            "model": {"method": "hf", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "h_mp2": {
        "program": "psi4",
        "specification": {
            "driver": "hessian",
            "model": {"method": "mp2", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "h_scf_atz": {
        "program": "psi4",
        "specification": {
            "driver": "hessian",
            "model": {"method": "hf", "basis": "aug-cc-pvtz"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "h_scf_adz": {
        "program": "psi4",
        "specification": {
            "driver": "hessian",
            "model": {"method": "hf", "basis": "aug-cc-pvdz"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
    "h_b3lyp": {
        "program": "psi4",
        "specification": {
            "driver": "hessian",
            "model": {"method": "b3lyp", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
}
