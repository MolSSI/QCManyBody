from pprint import pprint

from qcmanybody.manybody import Molecule
from qcmanybody.models import BsseEnum
from qcmanybody.qcengine import run_qcengine

mol = Molecule.from_data(
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

driver = "energy"
specs = {
    "scf": {
        "program": "psi4",
        "specification": {
            "driver": driver,
            "model": {"method": "hf", "basis": "sto-3g"},
            "keywords": {"cc_type": "df", "df_basis_mp2": "def2-qzvpp-ri"},
        },
    },
}

levels = {2: "scf"}

nbody_results = run_qcengine(mol, levels, specs, [BsseEnum.nocp], True)
pprint(nbody_results)
