# Core Interface

The core interface of QCManyBody is designed to allow for more flexibility in how the calculations are run.
The primary responsibilities of the core interface are:

1. Given a molecule and desired levels of MBE, return the fragments and levels that must be computed for each fragment
2. Given a dictionary of the results of those computations, analyze those results and calculate the desired manybody properties

Note that the user is expected to run the calculations themselves, and the core interface does not provide any tools for
running the calculations.


## Using the core interface

The core interface is accessed through the [`ManyBodyCalculator`][qcmanybody.manybody.ManyBodyCalculator]
class.

The first step is to create a molecule. This molecule is a QCElemental molecule object, and must contain fragments.

```python
import qcmanybody as qcmb
import qcelemental as qcel

# Create a molecule with 3 hydrogen atoms, each as its own fragment
mol = qcel.models.Molecule(symbols=['h', 'h', 'h'],
                           geometry=[[0,0,0],[0,0,2],[0,0,4]],
                           fragments=[[0], [1], [2]])
```

Next, create a `ManyBodyCalculator` object. This object is constructed using the molecule,
the desired BSSE correction, levels of MBE, and other options of the MBE.

```python
mbc = qcmb.ManyBodyCalculator(
    molecule=mol,
    bsse_type=['cp'],
    levels={1: 'mp2/aug-cc-pvtz', 2: 'b3lyp/def2-svp', 3: 'hf/sto-3g'},
    return_total_data=True,
    supersystem_ie_only=False
)
```

The `levels` option is a dictionary that specifies the n-body level as a key, then an arbitrary
string as the description of the calculation to be performed at the n-body level. This string is
termed the 'model chemistry' is completely arbitrary; it only has meaning to the user, and the user is expected to map these strings
to some meaningful definition of a calculation.

**Note**: The core interface is less flexible than the high-level interface when it comes to the `levels` option.
    In the core interface, all levels must be accounted for (that is, keys must go from 1 to the maximum
    nbody you would like to calculate). All levels must be present even if the model chemistry
    is the same for all levels.

For a complete discussion of the other options available in the `ManyBodyCalculator` object, see the
[keywords discussion](keywords.md)
the [`ManyBodyCalculator API documentation`][qcmanybody.manybody.ManyBodyCalculator].

The next step is to obtain the calculations to be run from the `ManyBodyCalculator` object.
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

Note that it is entirely up to the user to run the calculation somehow - this level of interface
does not provide any tools for running the calculations.

### Results dictionary

The data returned from the calculations is expected to be stored in a nested dictionary.
The level is the opaque label as given from the `QCManyBodyCalculator`.
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