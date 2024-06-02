# Keywords and options

## Required

### Molecule

The molecule used by QCManyBody is a [QCElemental Molecule object](https://molssi.github.io/QCElemental/model_molecule.html).
The only requirement for use in QCManyBody is that the molecule has multiple fragments. It is these fragments
that will be used in the many-body expansion.

Examples:

```python
from qcelemental.models import Molecule

# Molecule with 3 neon atoms, each as its own fragment
ne3 = Molecule(symbols=["ne", "ne", "ne"],
               geometry=[[0,0,0],[0,0,2],[0,0,4]],
               fragments=[[0], [1], [2]])

# Water tetramer
water4 = Molecule(symbols=['o', 'h', 'h', 'o', 'h', 'h', 'o', 'h', 'h', 'o', 'h', 'h'],
                  geometry=[[-2.8, -1.2, -0.2], [-1.1, -2.1, -0.0], [-3.8, -2.1,  1.1],
                            [-2.8, -1.2, -2.2], [-1.1, -2.1, -2.0], [-3.8, -2.1, -1.1],
                            [-2.8, -1.2, -4.2], [-1.1, -2.1, -4.0], [-3.8, -2.1, -3.1],
                            [-2.8, -1.2, -6.2], [-1.1, -2.1, -6.0], [-3.8, -2.1, -5.1]],
                          fragments=[[0,1,2], [3,4,5], [6,7,8], [9,10,11]])

# Water trimer, using the from_data function
# the -- is used to separate fragments
water3 = Molecule.from_data(
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

```

### bsse_type

The `bsse_type` parameter specifies the type of correction for basis set superposition error (BSSE). Multiple
types can be specified, in which case the [results](results.md) will include separate fields for each type of
correction.

Valid types are:

- `nocp` - No counterpoise or other corrections applied
- `cp` - Counterpoise correction
- `vmfc` - Valiron-Mayer function counterpoise correction

### levels and max_nbody

Dictionary of different levels of theory for different levels of expansion. The keys are integers or "supersystem",
and the values are arbitrary strings that represent the model chemistry or level of theory to use for that level.
This string is arbitrary and only has meaning to the user - the user is expected to map these strings to some
meaningful definition of a calculation.

If a `supersystem` key is present, all higher order n-body effects up to `max_nbody` will be computed.

In the [core interface](core-interface.md), all levels must be accounted for (that is, keys must go
from 1 to the maximum), and `max_nbody` cannot be specified. In the high-level interface, a computational model
fills in for any lower unlisted n-body levels.

In the high-level interface, if both levels and max_nbody are provided, they must be consistent.

Examples:

- `{2: 'ccsd(t)/cc-pvdz', 3: 'mp2'}`
- `max_nbody=3` and `levels={1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'}`


## Keywords and Options

### return_total_data

When set to true, the manybody calculation will return the total data (energy/gradient/hessian/property) of the system.
If not, the return will only contain interaction data.

Note that the calculation of counterpoise corrected total properties implies the calculation of the energies of monomers
in the monomer basis, hence specifying `return_total_data = True` may carry out more computations than.
For some properties such as gradients and hessians, `return_total_data = False` is rarely useful.

### supersystem_ie_only

Target the supersystem total/interaction energy (IE) data over the many-body expansion (MBE)
analysis, thereby omitting intermediate-body calculations. When false, each n-body level
in the MBE up through `max_nbody` will be computed. When true (only allowed for `max_nbody = nfragments`),
only compute enough for the overall interaction/total energy: max_nbody-body and 1-body.

When true, properties `INTERACTION {driver} THROUGH {max_nbody}-BODY` will always be available;
`TOTAL {driver} THROUGH {max_nbody}-BODY` will be available depending on `return_total_data`; and
`{max_nbody}-BODY CONTRIBUTION TO {driver}` won't be available (except for dimers).

This keyword produces no savings for a two-fragment molecule. But for the interaction energy of a three-fragment molecule, for example, 2-body
subsystems can be skipped with `supersystem_ie_only=True` Do not use with `vmfc` in `bsse_type`
as it cannot produce savings.
