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
    "default": {
        'driver': driver, 'model': {"method": "ccsd(t)", "basis": "aug-cc-pvdz"}
    },
}

mol = Molecule(symbols=["Ne", "Ne"], fragments=[[0], [1]], geometry=[0.0, 0.0, -2.834589188186742, 0.0, 0.0, 2.834589188186742])
levels = {1: "default", 2: "default"}

mc = ManybodyCalculator(mol, [BsseEnum.cp], levels, False)
print("COMPUTE MAP", mc.compute_map)

component_results = {'energy': {}, 'gradient': {}, 'hessian': {}}

for chem, label, imol in mc.iterate_molecules(True):
    filename = hashlib.sha256(f"{chem}_{mol.get_hash()}_{label}".encode()).hexdigest()
    filepath = f'{cache_dir}/{filename}.json.gz'

    if not os.path.exists(filepath):
        inp = AtomicInput(molecule=imol, **specs[chem])
        print(f"Computing {chem} | {label}")
        print(inp)
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
pprint(nbody_results)
