from __future__ import annotations

from typing import Dict, Union, Set

import numpy as np
from qcelemental.models import DriverEnum

from .models import FragBasIndex
from .utils import labeler


def _sum_cluster_ptype_data(
    ptype: DriverEnum,
    ptype_dict: Dict,
    compute_list: Set[FragBasIndex],
    fragment_slice_dict: Dict[int, slice],
    fragment_size_dict: Dict[int, int],
    mc_level_lbl: str,
    vmfc: bool = False,
    nb: int = 0,
) -> Union[float, np.ndarray]:
    """
    Sum arrays from n-body computations to obtain the BSSE corrected or uncorrected scalar or array.

    Parameters
    ----------
    ptype
        Hint to shape of array data to sum.
    ptype_dict
        Dictionary containing computed energy, gradient, or Hessian obtained from each subsystem computation
    compute_list
        A list of (frag, bas) tuples notating all the required computations.
    fragment_slice_dict
        Dictionary containing slices that index the gradient or Hessian matrix for each of the 1-indexed fragments.
        For He--HOOH--Me cluster, `{1: slice(0, 1, None), 2: slice(1, 5, None), 3: slice(5, 10, None)}`.
    fragment_size_dict
        Dictionary containing the number of atoms of each 1-indexed fragment.
        For He--HOOH--Me cluster, `{1: 1, 2: 4, 3: 5}`.
    vmfc
        Is it a VMFC calculation?
    nb
        n-body level; required for VMFC calculations.
    mc_level_lbl
        User label for what modelchem level results should be pulled out of *ptype_dict*.
        This is the 1-indexed counterpart to 0-indexed mc_level_idx.

    Returns
    -------
    ret
        Scalar or array containing the summed energy, gradient, or Hessian result.
        Formerly, passed in and modified in place and only called for g/h.

    """
    sign = 1
    nat = sum(fragment_size_dict.values())

    if ptype == DriverEnum.energy:
        ret = 0.0

        for frag, bas in compute_list:
            ene = ptype_dict[labeler(mc_level_lbl, frag, bas)]

            if vmfc:
                sign = (-1) ** (nb - len(frag))

            ret += sign * ene

        return ret

    elif ptype == DriverEnum.gradient:
        ret = np.zeros((nat, 3))

        for frag, bas in compute_list:
            grad = np.asarray(ptype_dict[labeler(mc_level_lbl, frag, bas)])

            if vmfc:
                sign = (-1) ** (nb - len(frag))

            start = 0
            for ifr in bas:
                end = start + fragment_size_dict[ifr]
                ret[fragment_slice_dict[ifr]] += sign * grad[start:end]
                start += fragment_size_dict[ifr]

        return ret

    elif ptype == DriverEnum.hessian:
        ret = np.zeros((nat * 3, nat * 3))

        for frag, bas in compute_list:
            hess = np.asarray(ptype_dict[labeler(mc_level_lbl, frag, bas)])

            if vmfc:
                sign = (-1) ** (nb - len(frag))

            # Build up start and end slices
            abs_start, rel_start = 0, 0
            abs_slices, rel_slices = [], []
            for ifr in bas:
                rel_end = rel_start + 3 * fragment_size_dict[ifr]
                rel_slices.append(slice(rel_start, rel_end))
                rel_start += 3 * fragment_size_dict[ifr]

                tmp_slice = fragment_slice_dict[ifr]
                abs_slices.append(slice(tmp_slice.start * 3, tmp_slice.stop * 3))

            for abs_sl1, rel_sl1 in zip(abs_slices, rel_slices):
                for abs_sl2, rel_sl2 in zip(abs_slices, rel_slices):
                    ret[abs_sl1, abs_sl2] += sign * hess[rel_sl1, rel_sl2]

        return ret

    else:
        raise KeyError(
            "ptype can only be energy, gradient, or hessian. How did you end up here?"
        )