import copy
import pprint
import re

import pytest
from qcelemental import constants
from qcelemental.models import Molecule

# v2: from qcelemental.models.procedures_manybody import AtomicSpecification, ManyBodyKeywords, ManyBodyInput
from qcelemental.testing import compare_recursive, compare_values

from qcmanybody.computer import ManyBodyComputer
from qcmanybody.models import AtomicSpecification, ManyBodyInput, ManyBodyKeywords, ManyBodyResultProperties
from qcmanybody.utils import translate_qcvariables

from .addons import using, uusing
from .test_computer_he4_singlelevel import sumdict as sumdict_single


@pytest.fixture(scope="function")
def mbe_data_multilevel_631g():
    # note that spherical/cartesian irrelevant for He & 6-31G, and fc/ae irrelevant for He
    c4_kwds = {}
    gms_kwds = {"basis__ngauss": 6, "ccinp__ncore": 0, "ccinp__iconv": 9, "scf__conv": 9}
    nwc_kwds = {"scf__thresh": 1.0e-8, "ccsd__thresh": 1.e-8}
    p4_kwds = {"scf_type": "pk", "mp2_type": "conv"}

    protocols = {"stdout": False}
    return {
        "specification": {
            "specification": {
                "c4-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "cfour",
                    "keywords": c4_kwds,
                    "protocols": protocols,
                },
                "c4-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "cfour",
                    "keywords": c4_kwds,
                    "protocols": protocols,
                },
                "c4-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "cfour",
                    "keywords": c4_kwds,
                    "protocols": protocols,
                },
                "gms-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "n31",
                    },
                    "driver": "energy",
                    "program": "gamess",
                    "keywords": gms_kwds,
                    "protocols": protocols,
                },
                "gms-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "n31",
                    },
                    "driver": "energy",
                    "program": "gamess",
                    "keywords": gms_kwds,
                    "protocols": protocols,
                },
                "gms-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "n31",
                    },
                    "driver": "energy",
                    "program": "gamess",
                    "keywords": gms_kwds,
                    "protocols": protocols,
                },
                "nwc-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "nwchem",
                    "keywords": nwc_kwds,
                    "protocols": protocols,
                },
                "nwc-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "nwchem",
                    "keywords": nwc_kwds,
                    "protocols": protocols,
                },
                "nwc-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "nwchem",
                    "keywords": nwc_kwds,
                    "protocols": protocols,
                },
                "p4-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "psi4",
                    "keywords": p4_kwds,
                    "protocols": protocols,
                },
                "p4-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "psi4",
                    "keywords": p4_kwds,
                    "protocols": protocols,
                },
                "p4-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "6-31g",
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


he4_refs_conv_multilevel_631g = {
    # 1: ccsd; 2,3: mp2; 4: hf, all 6-31G
    "121": {
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":       -11.480648555603,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":       -11.472000052247,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":       -11.472089645469,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":       -11.472068853166,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":   0.0,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   0.008648503357,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   0.008558910134,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   0.008579702437,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       0.008648503357,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":      -0.000089593222,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":       0.000020792303,

        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.480648555603,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.471058574581,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.471324608815,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.471272244751,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":     0.0,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.009589981022,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     0.009323946788,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.009376310852,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.009589981022,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        -0.000266034234,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         0.000052364064,

        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY":       -11.480648555603,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY":       -11.471163557706,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY":       -11.471402496106,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY":       -11.471350132042,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":   0.0,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   0.009484997897,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   0.009246059497,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   0.009298423561,
        "VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       0.009484997897,
        "VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":      -0.000238938400,
        "VMFC-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":       0.000052364064,
    },
    # 1,2: ccsd; 3,4: mp2, all 6-31G
    "22": {
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":       -11.480648555603,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":       -11.471764016410,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":       -11.471853609632,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":       -11.471834096023,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":   0.0,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   0.008884539193,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   0.008794945971,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   0.008814459580,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       0.008884539193,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":      -0.000089593222,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":       0.000019513609,

        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.480648555603,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.470705938773,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.470971973006,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.470913449084,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":     0.0,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.009942616831,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     0.009676582597,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.009735106519,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.009942616831,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        -0.000266034234,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         0.000058523922,

        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY":       -11.480648555603,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY":       -11.470821409457,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY":       -11.471060347857,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY":       -11.471001823935,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":   0.0,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   0.009827146147,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   0.009588207746,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   0.009646731668,
        "VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       0.009827146147,
        "VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":      -0.000238938400,
        "VMFC-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":       0.000058523922,
    },
    "ss22": {
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.480648555603,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.470705938773,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.470776016167,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY":     0.0,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.009942616831,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.009872539572,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.009942616831,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":        -0.000070079612,
    }
}

# only here for keys
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


sumstr = {
   "cp4b_tot": {
       "121": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+§B\s+2\s+-11.4710585\d+        0.0095899\d+        6.01\d+        0.0095899\d+        6.01\d+
^\s+§B\s+3\s+-11.4713246\d+        0.0093239\d+        5.85\d+       -0.0002660\d+       -0.16\d+
^\s+FULL/RTN\s+§C\s+4\s+-11.4712722\d+        0.0093763\d+        5.88\d+        0.0000523\d+        0.03\d+
""",
        "22": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+§A\s+2\s+-11.4707059\d+        0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
^\s+§B\s+3\s+-11.4709719\d+        0.0096765\d+        6.07\d+       -0.0002660\d+       -0.16\d+
^\s+FULL/RTN\s+§B\s+4\s+-11.4709134\d+        0.0097351\d+        6.10\d+        0.0000585\d+        0.03\d+
""",
        "ss22": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+§A\s+2\s+-11.4707059\d+        0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
^\s+SS\s+3\s+N/A\s+N/A\s+N/A\s+N/A\s+N/A\s*
^\s+SS/FULL/RTN\s+4      -11.4707760\d+        0.0098725\d+        6.19\d+       -0.0000700\d+       -0.04\d+
""",
   },
   "cp3b_tot": {
       "121": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+§B\s+2\s+-11.4710585\d+        0.0095899\d+        6.01\d+        0.0095899\d+        6.01\d+
^\s+RTN\s+§B\s+3\s+-11.4713246\d+        0.0093239\d+        5.85\d+       -0.0002660\d+       -0.16\d+
^\s+FULL\s+4\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
""",
        "22": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+§A\s+2\s+-11.4707059\d+        0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
^\s+RTN\s+§B\s+3\s+-11.4709719\d+        0.0096765\d+        6.07\d+       -0.0002660\d+       -0.16\d+
^\s+FULL\s+4\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
""",
   },
   "cp3b_ie": {
       "121": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+N/A\s+0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+§B\s+2\s+N/A\s+0.0095899\d+        6.01\d+        0.0095899\d+        6.01\d+
^\s+RTN\s+§B\s+3\s+N/A\s+0.0093239\d+        5.85\d+       -0.0002660\d+       -0.16\d+
^\s+FULL\s+4        N/A                   N/A                   N/A                   N/A                   N/A\s*
""",
        "22": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+N/A\s+0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+§A\s+2\s+N/A\s+0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
^\s+RTN\s+§B\s+3\s+N/A\s+0.0096765\d+        6.07\d+       -0.0002660\d+       -0.16\d+
^\s+FULL\s+4        N/A                   N/A                   N/A                   N/A                   N/A\s*
""",
   },
   "cp2b_tot": {
       "121": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+RTN\s+§B\s+2\s+-11.4710585\d+        0.0095899\d+        6.01\d+        0.0095899\d+        6.01\d+
^\s+3\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
^\s+FULL\s+4\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
""",
        "22": r"""
   ==> N-Body: Counterpoise Corrected \(CP\) energies <==

^\s+MC n-Body\s+Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
^\s+\[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
^\s+§A\s+1\s+-11.4806485\d+        0.0000000\d+        0.00\d+        0.0000000\d+        0.00\d+
^\s+RTN\s+§A\s+2\s+-11.4707059\d+        0.0099426\d+        6.23\d+        0.0099426\d+        6.23\d+
^\s+3\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
^\s+FULL\s+4\s+N/A                   N/A                   N/A                   N/A                   N/A\s*
""",
   },
}

sumdict_multi = {
    "4b_vmfc_rtd": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_vmfc": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",  # TODO remove?
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "3b_vmfc_rtd": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_vmfc": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY",  # TODO remove?
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "2b_vmfc_rtd": {
        "121": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        },
        "22": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        },
    },
    "2b_vmfc": {
        "121": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY",  # TODO remove?
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        },
        "22": {
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        },
    },
}
sumdict = copy.deepcopy(sumdict_single)
sumdict.update(sumdict_multi)


@pytest.fixture
def he_tetramer():
    a2 = 2 / constants.bohr2angstroms
    return Molecule(symbols=["He", "He", "He", "He"], fragments=[[0], [1], [2], [3]], geometry=[0, 0, 0, 0, 0, a2, 0, a2, 0, 0, a2, a2])


@uusing("qcengine")
@pytest.mark.parametrize("levels", [
    # pattern 121
    pytest.param({4: "c4-hf", 3: "c4-mp2", 1: "c4-ccsd"}, id="121-cfour_pure", marks=using("cfour")),
    pytest.param({4: "nwc-hf", 3: "nwc-mp2", 1: "nwc-ccsd"}, id="121-nwchem_pure", marks=using("nwchem")),
    pytest.param({4: "p4-hf", 3: "p4-mp2", 1: "p4-ccsd"}, id="121-psi4_pure", marks=using("psi4")),

    pytest.param({4: "p4-hf", 3: "c4-mp2", 1: "c4-ccsd"}, id="121-cfour_psi4", marks=[*using("cfour"), *using("psi4")]),
    pytest.param({4: "nwc-hf", 3: "nwc-mp2", 1: "p4-ccsd"}, id="121-nwchem_psi4", marks=[*using("nwchem"), *using("psi4")]),
    pytest.param({4: "c4-hf", 3: "nwc-mp2", 1: "p4-ccsd"}, id="121-cfour_nwchem_psi4", marks=[*using("cfour"), *using("nwchem"), *using("psi4")]),

    # pattern 22
    pytest.param({4: "c4-mp2", 2: "c4-ccsd"}, id="22-cfour_pure", marks=using("cfour")),
    pytest.param({4: "nwc-mp2", 2: "nwc-ccsd"}, id="22-nwchem_pure", marks=using("nwchem")),
    pytest.param({4: "p4-mp2", 2: "p4-ccsd"}, id="22-psi4_pure", marks=using("psi4")),
])
@pytest.mark.parametrize("mbe_keywords,anskey,bodykeys,outstrs,calcinfo_nmbe", [
#    pytest.param(
#        {"bsse_type": ["nocp", "cp", "vmfc"]},
#        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv],
#        {"121": 65,
#         "22": 99},  #
#        id="4b_all"),
    pytest.param(  # ODD entry b/c no vmfc ready for multilevel
        {"bsse_type": ["nocp", "cp"]},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if not k.startswith("VMFC-")],
        None,
        {"121": 61,  # cp(14md + 15lo) + nocp(14md + 15lo) + 4hi - 1lo 1234@1234
         "22": 49},  # cp(10hi + 15lo) + nocp(10hi + 15lo) - 1lo 1234@1234
        id="4b_nocpcp"),
#    pytest.param(
#        {"bsse_type": "nocp", "return_total_data": True, "supersystem_ie_only": True},
#        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
#        {"121": 5,
#         "22": 99},  #
#        id="4b_nocp_rtd_sio"),
#    pytest.param(
#        {"bsse_type": "nocp", "return_total_data": False, "supersystem_ie_only": True},
#        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
#        {"121": 5,
#         "22": 99},  #
#        id="4b_nocp_sio"),
#    pytest.param(
#        {"bsse_type": "cp", "return_total_data": True, "supersystem_ie_only": True},
#        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv if (k.startswith("CP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
#        {"121": 9,
#         "22": 99},  #
#        id="4b_cp_rtd_sio"),
#    pytest.param(
#        {"bsse_type": "cp", "return_total_data": False, "supersystem_ie_only": True},
#        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv if (k.startswith("CP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k) and "TOTAL ENERGY" not in k)],
#        {"121": 5,
#         "22": 99},  #
#        id="4b_cp_sio"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-"))],
        None,
        {"121": 33,  # 4hi + 14md + 15lo vs. 15 for single-level
         "22": 25},  # 10hi + 15lo
        id="4b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-"))],
        None,
        {"121": 33,  # could be 29 TODO,  # 14md + 15lo vs. 15 for single-level
         "22": 25},  # 10hi + 15lo
        id="4b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True},
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-"))],
        ["cp4b_tot"],
        {"121": 33,  # 4hi(nocp) + 14md + 15lo vs. 19 for single-level,
         "22": 29},  # 10hi + 15lo + 4hi(nocp)
        id="4b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "TOTAL ENERGY" not in k)],
        None,
        {"121": 29, # 14md + 15lo vs. 15 for single-level,
         "22": 25},  # 10hi + 15lo
        id="4b_cp"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": True},
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-"))],
        None,
        {"121": 65,  # was 93 in p4  # 4hi + 18+28md + 15lo vs. 65 for single-level
         "22": 65},  # was 83 in p4  # 4+18hi + 28+15lo
        id="4b_vmfc_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False},
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-"))],
        None,
        {"121": 65, # could be 61; was 93 in p4  # 18+28md + 15lo =61
         "22": 65}, # could be 61; was 83 in p4  # 18hi + 28+15lo = 61
        id="4b_vmfc"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 3},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k))],
        None,
        {"121": 18,  # 4hi + 14md vs. 14 for single-level
         "22": 24},  # 10hi + 14lo
        id="3b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 3},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and "4-BODY" not in k)],
        None,
        {"121": 18,  # 4hi + 14md vs. 14 for single-level
         "22": 24},  # 10hi + 14lo
        id="3b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 3},
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "4-BODY" not in k)],
        ["cp3b_tot"],
        {"121": 18,  # 4hi + 14md vs. 18 for single-level  # bugfix: was 28
         "22": 28},  # 10hi + 14lo + 4hi(nocp)
        id="3b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 3},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "4-BODY" not in k and "TOTAL ENERGY" not in k)],
        ["cp3b_ie"],
        {"121": 14,  # 14md vs. 14 for single-level
         "22": 24},  # 10hi + 14lo
        id="3b_cp"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": True, "max_nbody": 3},
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-") and ("4-BODY" not in k))],
        None,
        {"121": 50,  # 4hi + 18+28md vs. 14? for single-level
         "22": 50},  # was 68 in p4  # 4+18hi + 28lo
        id="3b_vmfc_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False, "max_nbody": 3},
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("VMFC-") and ("4-BODY" not in k))],
        None,
        {"121": 50, # could be 46  # 18+28md =46
         "22": 50}, # could be 46; was 68 in p4  # 18hi + 28lo = 46
        id="3b_vmfc"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 2},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        None,
        {"121": 14,  # 4hi + 10md vs. 10 for single-level
         "22": 10},  # 10hi
        id="2b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 2},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        None,
        {"121": 14,  # 4hi + 10md vs. 10 for single-level,
         "22": 10},  # 10hi
        id="2b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 2},
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        ["cp2b_tot"],
        {"121": 14,  # 4hi + 10md vs. 14 for single-level,
         "22": 14},  # 10hi + 4hi(nocp)
        id="2b_cp_rtd"),
    pytest.param(
        {"bsse_type": "ssfc", "return_total_data": False, "max_nbody": 2},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("4-BODY" not in k) and ("3-BODY" not in k) and "TOTAL ENERGY" not in k)],
        None,
        {"121": 10,  # 10md vs. 10 for single-level,
         "22": 10},  # 10hi
        id="2b_cp"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": True, "max_nbody": 2},
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        # "22" at 2-body is effectively single-level so nocp enabled ...
        {"121": [k for k in he4_refs_conv if (k.startswith("VMFC-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
         "22": [k for k in he4_refs_conv if ((k.startswith("VMFC-") or k.startswith("NOCP-")) and ("4-BODY" not in k) and ("3-BODY" not in k))]},
         None,
        {"121": 22,  # 4hi + 18+28md + 15lo vs. 65 for single-level
         "22": 22},  # 4+18hi + 28+15lo
        id="2b_vmfc_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False, "max_nbody": 2},
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        {"121": [k for k in he4_refs_conv if (k.startswith("VMFC-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
         "22": [k for k in he4_refs_conv if ((k.startswith("VMFC-") or k.startswith("NOCP-")) and ("4-BODY" not in k) and ("3-BODY" not in k))]},
         None,
        {"121": 22, # TODO could be 18  # 0+18hi = 18
         "22": 22}, # TODO could be 18  # 18hi = 18
        id="2b_vmfc"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 1},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("1-BODY" in k))],
        None,
        {"121": 4,  # 4hi
         "22": 4},
        id="1b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 1},
       "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
       [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("1-BODY" in k))],
       None,
        {"121": 4,  # 4hi
         "22": 4},
        id="1b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 1},
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("1-BODY" in k))],
        None,
        {"121": 4,  # 4hi
         "22": 4},
        id="1b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 1},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("1-BODY" in k) and "TOTAL ENERGY" not in k)],
        None,
        {"121": 0,
         "22": 0},
        id="1b_cp"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": True, "max_nbody": 1},
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        # max_nbody=1 is effectively single-modelchem so NOCP available for free
        [k for k in he4_refs_conv if ((k.startswith("VMFC-") or k.startswith("NOCP-")) and ("1-BODY" in k))],
        None,
        {"121": 4,  # 4hi
         "22": 4},
        id="1b_vmfc_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False, "max_nbody": 1},
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if ((k.startswith("VMFC-") or k.startswith("NOCP-")) and ("1-BODY" in k))],
        None,
        {"121": 4, # TODO could be 0
         "22": 4}, # TODO could be 0
        id="1b_vmfc"),
])
def test_nbody_he4_multi(levels, mbe_keywords, anskey, bodykeys, outstrs, calcinfo_nmbe, he_tetramer, request, mbe_data_multilevel_631g):
    _inner = request.node.name.split("[")[1].split("]")[0]
    kwdsln, pattern, progln = _inner.split("-")

    levels = copy.deepcopy(levels)
    if pattern == "121":
        if mbe_keywords.get("max_nbody", 4) == 3:
            del levels[4]  # max_nbody and levels silently contradict w/o this
        elif mbe_keywords.get("max_nbody", 4) == 2:
            levels = {2: levels[3], 1: levels[1]}
        elif mbe_keywords.get("max_nbody", 4) == 1:
            del levels[4]
            del levels[3]
    elif pattern == "22":
        if mbe_keywords.get("max_nbody", 4) == 3:
            levels = {3: levels[4], 2: levels[2]}
        if mbe_keywords.get("max_nbody", 4) == 2:
            del levels[4]
        if mbe_keywords.get("max_nbody", 4) == 1:
            levels = {1: levels[2]}

    mbe_keywords = ManyBodyKeywords(levels=levels, **mbe_keywords)
    mbe_data_multilevel_631g["molecule"] = he_tetramer
    mbe_data_multilevel_631g["specification"]["keywords"] = mbe_keywords
    mbe_model = ManyBodyInput(**mbe_data_multilevel_631g)

    # qcng: ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    ret = ManyBodyComputer.from_manybodyinput(mbe_model)
    print(f"MMMMMMM {request.node.name}")
    pprint.pprint(ret.dict(), width=200)

    # don't want QCVariables stashed in extras, but prepare the qcvars translation, and check it
    assert ret.extras == {}, f"[w] extras wrongly present: {ret.extras.keys()}"
    qcvars = translate_qcvariables(ret.properties.dict())

    skprop = ManyBodyResultProperties.to_qcvariables(reverse=True)

    refs = he4_refs_conv_multilevel_631g[pattern]
    ans = refs[anskey]
    ref_nmbe = calcinfo_nmbe[pattern]
    ref_bodykeys = bodykeys[pattern] if pattern in bodykeys else bodykeys
    ref_sumdict = sumdict[kwdsln][pattern] if pattern in sumdict[kwdsln] else sumdict[kwdsln]
    atol = 2.5e-8

    for qcv, ref in refs.items():
        skp = skprop[qcv]
        if qcv in ref_bodykeys:
            assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[b] skprop {skp}")
        else:
            assert qcv not in qcvars, f"[z] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv in sumdict["4b_all"]:
        skp = skprop[qcv]
        if qcv in ref_sumdict:
            ref = refs[ref_sumdict[qcv]]
            assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[c] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[d] skprop {skp}")
        else:
            assert qcv not in qcvars, f"[y] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop[qcv]
        assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=atol, label=f"[g] ret")

    assert ret.properties.calcinfo_nmbe == ref_nmbe, f"{ret.properties.calcinfo_nmbe=} != {ref_nmbe}"

    if outstrs:
        for stdoutkey in outstrs:
            assert re.search(sumstr[stdoutkey][pattern], ret.stdout, re.MULTILINE), f"[j] N-Body pattern not found: {sumstr[stdoutkey][pattern]}"


@uusing("qcengine")
@pytest.mark.parametrize("levels", [
    # pattern ss121
    #pytest.param({4: "c4-hf", 3: "c4-mp2", 1: "c4-ccsd"}, id="121-cfour_pure", marks=using("cfour")),
    #pytest.param({4: "nwc-hf", 3: "nwc-mp2", 1: "nwc-ccsd"}, id="121-nwchem_pure", marks=using("nwchem")),
    #pytest.param({"supersystem": "p4-hf", 3: "p4-mp2", 1: "p4-ccsd"}, id="ss121-psi4_pure", marks=using("psi4")),

    # pattern ss22
    #pytest.param({4: "c4-mp2", 2: "c4-ccsd"}, id="22-cfour_pure", marks=using("cfour")),
    #pytest.param({4: "nwc-mp2", 2: "nwc-ccsd"}, id="22-nwchem_pure", marks=using("nwchem")),
    pytest.param({"supersystem": "p4-mp2", 2: "p4-ccsd"}, id="ss22-psi4_pure", marks=using("psi4")),
])
@pytest.mark.parametrize("mbe_keywords,anskey,bodykeys,outstrs,calcinfo_nmbe", [
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True},
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-"))],
        ["cp4b_tot"],
        {#"ss121": 0,
         "ss22": 25},  # cp(10hi) + nocp(4hi + 11lo)
        id="4b_cp_rtd"),
])
def test_nbody_he4_supersys(levels, mbe_keywords, anskey, bodykeys, outstrs, calcinfo_nmbe, he_tetramer, request, mbe_data_multilevel_631g):
    _inner = request.node.name.split("[")[1].split("]")[0]
    kwdsln, pattern, progln = _inner.split("-")

    levels = copy.deepcopy(levels)
    if pattern == "121":
        if mbe_keywords.get("max_nbody", 4) == 3:
            del levels[4]  # max_nbody and levels silently contradict w/o this
        elif mbe_keywords.get("max_nbody", 4) == 2:
            levels = {2: levels[3], 1: levels[1]}
        elif mbe_keywords.get("max_nbody", 4) == 1:
            del levels[4]
            del levels[3]
    elif pattern == "22":
        if mbe_keywords.get("max_nbody", 4) == 3:
            levels = {3: levels[4], 2: levels[2]}
        if mbe_keywords.get("max_nbody", 4) == 2:
            del levels[4]
        if mbe_keywords.get("max_nbody", 4) == 1:
            levels = {1: levels[2]}

    mbe_keywords = ManyBodyKeywords(levels=levels, **mbe_keywords)
    mbe_data_multilevel_631g["molecule"] = he_tetramer
    mbe_data_multilevel_631g["specification"]["keywords"] = mbe_keywords
    mbe_model = ManyBodyInput(**mbe_data_multilevel_631g)

    # qcng: ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    ret = ManyBodyComputer.from_manybodyinput(mbe_model)
    print(f"MMMMMMM {request.node.name}")
    pprint.pprint(ret.dict(), width=200)

    # don't want QCVariables stashed in extras, but prepare the qcvars translation, and check it
    assert ret.extras == {}, f"[w] extras wrongly present: {ret.extras.keys()}"
    qcvars = translate_qcvariables(ret.properties.dict())

    skprop = ManyBodyResultProperties.to_qcvariables(reverse=True)

    refs = he4_refs_conv_multilevel_631g[pattern]
    ans = refs[anskey]
    ref_nmbe = calcinfo_nmbe[pattern]
    ref_bodykeys = bodykeys[pattern] if pattern in bodykeys else bodykeys
    ref_sumdict = sumdict[kwdsln][pattern] if pattern in sumdict[kwdsln] else sumdict[kwdsln]
    atol = 2.5e-8

    for qcv, ref in refs.items():
        skp = skprop[qcv]
        if qcv in ref_bodykeys:
            assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[b] skprop {skp}")
        else:
            assert qcv not in qcvars, f"[z] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv in sumdict["4b_all"]:
        skp = skprop[qcv]
        if qcv in ref_sumdict:
            ref = refs[ref_sumdict[qcv]]
            assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[c] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[d] skprop {skp}")
        else:
            assert qcv not in qcvars, f"[y] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop[qcv]
        assert compare_values(ref, qcvars[qcv], atol=atol, label=f"[e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=atol, label=f"[g] ret")

    assert ret.properties.calcinfo_nmbe == ref_nmbe, f"{ret.properties.calcinfo_nmbe=} != {ref_nmbe}"

    if outstrs:
        for stdoutkey in outstrs:
            print(stdoutkey)
            print(sumstr[stdoutkey][pattern])
            assert re.search(sumstr[stdoutkey][pattern], ret.stdout, re.MULTILINE), f"[j] N-Body pattern not found: {sumstr[stdoutkey][pattern]}"
