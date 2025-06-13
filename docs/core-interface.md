# Core Interface

The core interface of QCManyBody is designed to allow for more flexibility in how the calculations are run.
The primary responsibilities of the core interface are:

1. Given a molecule and desired levels of MBE, return the fragments and levels that must be computed for each fragment
2. Given a dictionary of the results of those computations, analyze those results and calculate the desired many-body properties

Note that the user is expected to run the calculations themselves, and the core interface does not provide any tools for
running the calculations.


## Using the core interface

The core interface is accessed through the [`ManyBodyCore`][qcmanybody.core.ManyBodyCore]
class.

The first step is to create a molecule. This molecule is a
[QCElemental molecule object](https://molssi.github.io/QCElemental/model_molecule.html), and must contain fragments.
(see also: [moleule input](keywords.md#molecule))

```python
from qcelemental.models import Molecule

# Create a molecule with 3 neon atoms, each as its own fragment

mol = Molecule(symbols=["ne", "ne", "ne"],
               geometry=[[0,0,0], [0,0,2], [0,0,4]],
               fragments=[[0], [1], [2]])
```

Next, create a `ManyBodyCore` object. This object is constructed using the molecule,
the desired BSSE correction, levels of MBE, and other options of the MBE.

```python
from qcmanybody import ManyBodyCore

mbc = ManyBodyCore(molecule=mol,
                   bsse_type=["cp"],
                   levels={1: "ccsd/cc-pvtz",
                           2: "mp2/cc-pvdz",
                           3: "mp2/cc-pvdz"},
                   return_total_data=True,
                   supersystem_ie_only=False,
                   embedding_charges=None)
```

The `levels` option is a dictionary that specifies the n-body level as a key, then an arbitrary
string as the description of the calculation to be performed at the n-body level. This string is
termed the 'model chemistry' is completely arbitrary; it only has meaning to the user, and the user is expected to
map these strings to some meaningful definition of a calculation.

**Note:** The core interface is less flexible than the high-level interface when it comes to the `levels` option.
    In the core interface, all levels must be accounted for (that is, keys must go from 1 to the maximum
    nbody you would like to calculate). All levels must be present even if the model chemistry
    is the same for all levels.

For a complete discussion of the other options available in the `ManyBodyCore` object, see the
[keywords discussion](keywords.md)
the [`ManyBodyCore API documentation`][qcmanybody.core.ManyBodyCore].

The next step is to obtain the calculations to be run from the `ManyBodyCore` object.
This is done with a python generator function `iterate_molecules` that returns
a tuple. This tuple contains

1. The string describing the calculation to be run (the model chemistry string, as defined in the `levels` dictionary)
2. A label for the calculation. This label is opaque to the user but used to identify the calculation when analyzing the results.
3. A `Molecule` object that contains the cluster of fragments to be computed.

```python
calculation_results = {}
for model_chemistry, label, mol_cluster in mbc.iterate_molecules():
    calculation_results[label] = run_calculation(mol_cluster, model_chemistry)
```

Note that it is entirely up to the user to run the calculation somehow (`run_calculation`
above is not a provided function) - this level of interface
does not provide any tools for running the calculations.
See [test_example.py](https://github.com/MolSSI/QCManyBody/blob/main/qcmanybody/tests/test_examples.py)
for a worked example or the below as a quick example of what `run_calculation` would look like using QCEngine.

```python
import qcelemental as qcel
import qcengine as qcng
import qcmanybody as qcmb

p4_string = """
He 0 0 0
--
Ne 1 0 0
"""
bsse_type = "cp"
method = "mp2"
basis = "def2-svp"
program = "psi4"

mol = qcel.models.Molecule.from_data(p4_string)
mbc = qcmb.ManyBodyCore(
    molecule=mol,
    bsse_type=[bsse_type],
    levels={
        1: f"{method}/{basis}",
        2: f"{method}/{basis}",
    },
    return_total_data=True,
    supersystem_ie_only=False,
    embedding_charges=None,
)
calculation_results = {}
for chem, label, imol in mbc.iterate_molecules():
    mc, real, bas = qcmb.utils.delabeler(label)
    atin = qcel.models.AtomicInput(
        molecule=imol, driver="energy", model={"method": method, "basis": basis}
    )
    atres = qcng.compute(atin, program, raise_error=True)
    calculation_results[label] = {
        "energy": atres.properties.return_energy,
    }

final_results = mbc.analyze(calculation_results)
IE = final_results["results"][f"{bsse_type.lower()}_corrected_interaction_energy"]
print(IE)
#> 0.409787
```


### Results dictionary

The data returned from the calculations is expected to be stored in a nested dictionary.
The level is the opaque label as given from the `QCManyBodyCore`.
The second level is the name of the property.

```python
calculation_results = {
    'label1': {
        'energy': -1.0,
        'gradient': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    },
    'label2': {
        'energy': -2.0,
        'gradient': [[4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
    },
    'label3': {
        'energy': -3.0,
        'gradient': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
    },
}
```

### Analysis

This result dictionary is all that is needed to perform the final analysis and calculation
of the MBE properties.

```python
final_results = mbc.analyze(calculation_results)
```

For a discussion about what the results contain, see the [results documentation](results.md).
