from __future__ import annotations

import json
from typing import Tuple, Dict, Union, Iterable

import numpy as np
from qcelemental import constants


def zeros_like(x: Union[float, np.ndarray]) -> Union[int, float, np.ndarray]:
    if isinstance(x, float):
        return 0.0
    else:
        return np.zeros_like(x)


def copy_value(x: Union[float, np.ndarray]) -> Union[int, float, np.ndarray]:
    if isinstance(x, float):
        return x
    else:
        return np.copy(x)


def all_same_shape(it: Iterable[Union[float, np.ndarray]]) -> bool:
    """Check if all elements of an iterable have the same shape."""

    it = iter(it)
    try:
        first = next(it)
    except StopIteration:
        return True
    if isinstance(first, float):
        return all(isinstance(x, float) for x in it)
    elif isinstance(first, np.ndarray):
        return all(x.shape == first.shape for x in it)
    else:
        raise TypeError(f"Expected float or np.ndarray, got {type(first)}")


def expand_gradient(
    grad: np.ndarray, bas: Tuple[int, ...], fragment_size_dict: Dict[int, int], fragment_slice_dict: Dict[int, slice]
) -> np.ndarray:
    """
    Expands a gradient calculated for a cluster to the full system
    """

    nat = sum(fragment_size_dict.values())
    ret = np.zeros((nat, 3))
    start = 0
    for ifr in bas:
        end = start + fragment_size_dict[ifr]
        ret[fragment_slice_dict[ifr]] = grad[start:end]
        start += fragment_size_dict[ifr]

    return ret


def expand_hessian(
    hess: np.ndarray, bas: Tuple[int, ...], fragment_size_dict: Dict[int, int], fragment_slice_dict: Dict[int, slice]
) -> np.ndarray:
    """
    Expands a hessian calculated for a cluster to the full system
    """

    nat = sum(fragment_size_dict.values())
    ret = np.zeros((nat * 3, nat * 3))

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
            ret[abs_sl1, abs_sl2] = hess[rel_sl1, rel_sl2]

    return ret


def labeler(mc_level_lbl: str, frag: Tuple[int, ...], bas: Tuple[int, ...]) -> str:
    return json.dumps([str(mc_level_lbl), frag, bas])


def delabeler(item: str) -> Tuple[str, Tuple[int, ...], Tuple[int, ...]]:
    """Transform labels like string "1_((2,), (1, 2))" into tuple (1, (2,), (1, 2))."""

    mc, frag, bas = json.loads(item)
    return str(mc), frag, bas


def print_nbody_energy(
    energy_body_dict: Dict[int, float],
    header: str,
    nfragments: int,
    embedding: bool = False,
):
    """Format output string for user for a single bsse_type. Prints to output and logger.
    Called repeatedly by assemble_nbody_component."""

    info = f"""\n   ==> N-Body: {header} energies <==\n\n"""
    info += f"""  {"n-Body":>12}     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy\n"""
    info += f"""                   [Eh]                    [Eh]                  [kcal/mol]            [Eh]                  [kcal/mol]\n"""
    previous_e = energy_body_dict[1]
    tot_e = previous_e != 0.0
    nbody_range = list(energy_body_dict)
    nbody_range.sort()
    for nb in range(1, nfragments + 1):
        lbl = []
        if nb == nfragments:
            lbl.append("FULL")
        if nb == max(nbody_range):
            lbl.append("RTN")
        lbl = "/".join(lbl)

        if nb in nbody_range:
            delta_e = energy_body_dict[nb] - previous_e
            delta_e_kcal = delta_e * constants.hartree2kcalmol
            if embedding:
                int_e = np.nan
                int_e_kcal = np.nan
            else:
                int_e = energy_body_dict[nb] - energy_body_dict[1]
                int_e_kcal = int_e * constants.hartree2kcalmol
            if tot_e:
                info += f"""  {lbl:>8} {nb:3}  {energy_body_dict[nb]:20.12f}  {int_e:20.12f}  {int_e_kcal:20.12f}  {delta_e:20.12f}  {delta_e_kcal:20.12f}\n"""
            else:
                info += f"""  {lbl:>8} {nb:3}  {"N/A":20}  {int_e:20.12f}  {int_e_kcal:20.12f}  {delta_e:20.12f}  {delta_e_kcal:20.12f}\n"""
            previous_e = energy_body_dict[nb]
        else:
            info += f"""  {lbl:>8} {nb:3}        {"N/A":20}  {"N/A":20}  {"N/A":20}  {"N/A":20}  {"N/A":20}\n"""

    info += "\n"
    print(info)


def collect_vars(bsse, body_dict, max_nbody: int, embedding: bool = False):
    previous_e = body_dict[1]
    tot_e = previous_e != 0.0
    nbody_range = list(body_dict)
    nbody_range.sort()
    res = {}
    if embedding:
        return res

    if tot_e:
        res[f"{bsse}-CORRECTED TOTAL ENERGY"] = body_dict[max_nbody]
    res[f"{bsse}-CORRECTED INTERACTION ENERGY"] = body_dict[max_nbody] - body_dict[1]

    for nb in range(2, max(nbody_range) + 1):
        res[f"{bsse}-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY"] = body_dict[nb] - body_dict[1]
        res[f"{bsse}-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY"] = body_dict[nb] - body_dict[nb - 1]
    if tot_e:
        for nb in nbody_range:
            res[f"{bsse}-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY"] = body_dict[nb]

    return res
