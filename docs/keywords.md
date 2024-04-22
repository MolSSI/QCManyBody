# Keywords and options

## Required

### Molecule

### bsse_type

### levels and max_nbody

## Keywords and Options

Both the high-level interface and core interface share the same terminology with respect to options


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
