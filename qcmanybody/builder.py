from __future__ import annotations

import itertools
from typing import Iterable, Union, Literal, Optional, Dict, Set

from qcmanybody.models import BsseEnum, FragBasIndex


def build_nbody_compute_list(
    bsse_type: Iterable[BsseEnum],
    nfragments: int,
    nbodies: Iterable[Union[int, Literal["supersystem"]]],
    return_total_data: bool,
    supersystem_max_nbody: Optional[int] = None,
) -> Dict[str, Dict[int, Set[FragBasIndex]]]:
    """Generates lists of N-Body computations needed for requested BSSE treatments.

    Parameters
    ----------
    bsse_type
        Requested BSSE treatments.
    nfragments
        Number of distinct fragments comprising the full molecular supersystem.
    nbodies
        List of n-body levels (e.g., `[2]` or `[1, 2]` or `["supersystem"]`) for which to generate tasks.
        Note the natural 1-indexing, so `[1]` covers one-body contributions.
    return_total_data
        Whether the total data (True; energy/gradient/Hessian) of the molecular system has been requested,
        as opposed to interaction data (False).

    Returns
    -------
    compute_dict
        Dictionary containing subdicts enumerating compute lists for each possible BSSE treatment.
        Subdict keys are n-body levels and values are sets of all the `mc_(frag, bas)` indices
        needed to compute that n-body level. A given index can appear multiple times within a
        subdict and among subdicts.

            compute_dict["cp"] = {
                1: set(),
                2: {((1,), (1, 2)),
                    ((2,), (1, 2)),
                    ((1, 2), (1, 2))}
            }

        Subdicts below are always returned. Any may be empty if not requested through *bsse_type*.

        * ``'all'`` |w---w| full list of computations required
        * ``'cp'`` |w---w| list of computations required for CP procedure
        * ``'nocp'`` |w---w| list of computations required for non-CP procedure
        * ``'vmfc_compute'`` |w---w| list of computations required for VMFC procedure
        * ``'vmfc_levels'`` |w---w| list of levels required for VMFC procedure

    """

    include_supersystem = False
    if "supersystem" in nbodies:
        if supersystem_max_nbody is None:
            raise ValueError("supersystem_max_nbody must be provided if 'supersystem' contains nbodies")

        include_supersystem = True
        nbodies = list(nbodies)
        nbodies.remove("supersystem")

    # What levels do we need?
    fragment_range = range(1, nfragments + 1)

    # Need nbodies and all lower-body in full basis
    cp_compute_list = {x: set() for x in nbodies}
    nocp_compute_list = {x: set() for x in nbodies}
    vmfc_compute_list = {x: set() for x in nbodies}
    vmfc_level_list = {x: set() for x in nbodies}  # Need to sum something slightly different

    # Verify proper passing of bsse_type. already validated in Computer
    bsse_type_remainder = set(bsse_type) - {e.value for e in BsseEnum}
    if bsse_type_remainder:
        raise RuntimeError(f"Unrecognized BSSE type(s): {bsse_type_remainder}")

    # Build up compute sets
    if "cp" in bsse_type:
        # Everything is in full n-mer basis
        basis_tuple = tuple(fragment_range)

        for nb in nbodies:
            if nb > 1:
                for sublevel in range(1, nb + 1):
                    for x in itertools.combinations(fragment_range, sublevel):
                        cp_compute_list[nb].add((x, basis_tuple))

    if "nocp" in bsse_type:
        # Everything in monomer basis
        for nb in nbodies:
            for sublevel in range(1, nb + 1):
                for x in itertools.combinations(fragment_range, sublevel):
                    nocp_compute_list[nb].add((x, x))

    if "vmfc" in bsse_type:
        # Like a CP for all combinations of pairs or greater
        for nb in nbodies:
            for cp_combos in itertools.combinations(fragment_range, nb):
                basis_tuple = tuple(cp_combos)
                for interior_nbody in range(1, nb + 1):
                    for x in itertools.combinations(cp_combos, interior_nbody):
                        combo_tuple = (x, basis_tuple)
                        vmfc_compute_list[nb].add(combo_tuple)
                        vmfc_level_list[len(basis_tuple)].add(combo_tuple)

    if return_total_data and 1 in nbodies:
        # Monomers in monomer basis
        nocp_compute_list.setdefault(1, set())
        for frag in range(1, nfragments + 1):
            nocp_compute_list[1].add(((frag,), (frag,)))

    if include_supersystem:
        # Add supersystem info to the compute list (nocp only)
        for nb in range(1, supersystem_max_nbody + 1):
            cp_compute_list.setdefault(nb, set())
            nocp_compute_list.setdefault(nb, set())
            vmfc_compute_list.setdefault(nb, set())
            for sublevel in range(1, nb + 1):
                for x in itertools.combinations(fragment_range, sublevel):
                    nocp_compute_list[nb].add((x, x))

        # Add the total supersystem (nfragments@nfragments)
        nocp_compute_list.setdefault(nfragments, set())
        nocp_compute_list[nfragments].add((tuple(fragment_range), tuple(fragment_range)))

    # Build a comprehensive compute range
    # * do not use list length to count number of {nb}-body computations
    compute_list = {x: set() for x in nbodies}
    for nb in nbodies:
        compute_list[nb] |= cp_compute_list[nb]
        compute_list[nb] |= nocp_compute_list[nb]
        compute_list[nb] |= vmfc_compute_list[nb]

    if include_supersystem:
        for nb, lst in nocp_compute_list.items():
            compute_list.setdefault(nb, set())
            compute_list[nb] |= lst

    # Rearrange compute_list from key nb having values to compute all of that nb
    #   to key nb including values of that nb. Use for counting.
    compute_list_count = {x: set() for x in nbodies}
    for nb in nbodies:
        for nbset in compute_list.values():
            for item in nbset:
                if len(item[0]) == nb:
                    compute_list_count[nb].add(item)

    compute_dict = {
        "all": compute_list,
        "cp": cp_compute_list,
        "nocp": nocp_compute_list,
        "vmfc_compute": vmfc_compute_list,
        "vmfc_levels": vmfc_level_list,
    }
    return compute_dict
