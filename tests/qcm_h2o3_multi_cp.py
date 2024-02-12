from pprint import pprint

import os
import hashlib
import json
import gzip
import qcengine as qcng
from qcelemental.models import AtomicInput, AtomicResult

from qcmanybody.manybody import Molecule, ManybodyCalculator
from qcmanybody.models import BsseEnum

cache_dir = "/tmp/qcmanybodycache"
os.makedirs(cache_dir, exist_ok=True)

driver = "energy"
specs = {
    "scf": {
        'driver': driver,
        'model': {"method": "hf", "basis": "sto-3g"},
        'keywords': {'cc_type': 'df', 'df_basis_mp2': "def2-qzvpp-ri"},
    },
    "mp2": {
        'driver': driver,
        'model': {"method": "mp2", "basis": "sto-3g"},
        'keywords': {'cc_type': 'df', 'df_basis_mp2': "def2-qzvpp-ri"},
    },
}

mol = Molecule.from_data("""
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
""")

levels = {2: "scf", 1: "mp2"}

mc = ManybodyCalculator(mol, [BsseEnum.cp], levels, True)

component_results = {'energy': {}, 'gradient': {}, 'hessian': {}}

for chem, label, imol in mc.iterate_molecules(True):
    filename = hashlib.sha256(f"{chem}_{mol.get_hash()}_{label}".encode()).hexdigest()
    filepath = f'{cache_dir}/{filename}.json.gz'

    if not os.path.exists(filepath):
        inp = AtomicInput(molecule=imol, **specs[chem])
        print(f"Computing {chem} | {label}")
        result = qcng.compute(inp, "psi4")

        if not result.success:
            print(result.error.error_message)
            raise RuntimeError("Not successful")

        with gzip.open(filepath, "wt") as f:
            json.dump(result.dict(encoding='json'), f)

    with gzip.open(filepath, "rt") as f:
        result = AtomicResult(**json.load(f))

    component_results['energy'][label] = result.properties.return_energy
    component_results['gradient'][label] = result.properties.return_gradient
    component_results['hessian'][label] = result.properties.return_hessian


nbody_results = mc.analyze(driver, component_results['energy'])
pprint(component_results)
pprint(nbody_results)
