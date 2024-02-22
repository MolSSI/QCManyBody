from pprint import pprint

from common import mol_h2o_3, specifications
from qcmanybody.models import BsseEnum
from qcmanybody.qcengine import run_qcengine

#levels = {1: "e_scf", 2: "e_scf", 3: "e_scf"}
#levels = {1: "g_scf"}#, 2: "g_scf", 3: "g_scf"}
#levels = {1: "h_scf", 2: "h_scf", 3: "h_scf"}

levels = {1: "e_mp2", 2: "e_b3lyp", "supersystem": "e_scf"}

nbody_results = run_qcengine(mol_h2o_3, levels, specifications,
                             [BsseEnum.cp, BsseEnum.nocp, BsseEnum.vmfc], True)
pprint(nbody_results)