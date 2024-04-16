from __future__ import annotations

import json
from typing import Dict, Mapping, Tuple, Union, Iterable

import numpy as np
from qcelemental import constants


def find_shape(x: Union[float, np.ndarray]) -> Tuple[int, ...]:
    if isinstance(x, float):
        return (1,)
    else:
        return x.shape


def shaped_zero(shape: Tuple[int, ...]) -> Union[float, np.ndarray]:
    if shape == (1,):
        return 0.0
    else:
        return np.zeros(shape)


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


def resize_gradient(
    grad: np.ndarray,
    bas: Tuple[int, ...], 
    fragment_size_dict: Dict[int, int], 
    fragment_slice_dict: Dict[int, slice], 
    *, 
    reverse: bool = False,
) -> np.ndarray:
    """Pads or extracts a gradient array between subsystem and full supersystem sizes.

    Parameters
    ----------
    grad
        Gradient matrix of natural size for *bas*, (3 * <nat in bas>, 3).
        If `reverse=True`, gradient matrix of supersystem size, (3 * <nat of all fragments>, 3).
    bas
        1-indexed fragments active in *grad*.
        If `reverse=True`, 1-indexed fragments to be extracted from *grad*.
    fragment_size_dict
        Dictionary containing the number of atoms of each 1-indexed fragment.
        For He--HOOH--Me cluster, `{1: 1, 2: 4, 3: 5}`.
    fragment_slice_dict
        Dictionary containing slices that index the gradient matrix for each of the 1-indexed fragments.
        For He--HOOH--Me cluster, `{1: slice(0, 1), 2: slice(1, 5), 3: slice(5, 10)}`.
    reverse
        If True, contract *grad* from supersystem size and shape that which is natural for *bas*.

    Returns
    -------
        Gradient array padded with zeros to supersystem size, (3 * <nat of supersystem>, 3).
        If reverse=True, gradient array extracted to natural size, (3 * <nat in bas>, 3).

    """
    if reverse:
        nat = sum(fragment_size_dict[ifr] for ifr in bas)
    else:
        nat = sum(fragment_size_dict.values())
    ret = np.zeros((nat, 3))

    start = 0
    for ifr in bas:
        end = start + fragment_size_dict[ifr]
        if reverse:
            ret[start:end] = grad[fragment_slice_dict[ifr]]
        else:
            ret[fragment_slice_dict[ifr]] = grad[start:end]
        start += fragment_size_dict[ifr]

    return ret


def resize_hessian(
    hess: np.ndarray,
    bas: Tuple[int, ...],
    fragment_size_dict: Dict[int, int],
    fragment_slice_dict: Dict[int, slice],
    *,
    reverse: bool = False,
) -> np.ndarray:
    """Pads or extracts a Hessian array between subsystem and full supersystem sizes.

    Parameters
    ----------
    grad
        Hessian matrix of natural size for *bas*, (3 * <nat in bas>, 3 * <nat in bas>).
        If `reverse=True`, Hessian matrix of supersystem size, (3 * <nat of all fragments>, 3 * <nat of all fragments>).
    bas
        1-indexed fragments active in *hess*.
        If `reverse=True`, 1-indexed fragments to be extracted from *hess*.
    fragment_size_dict
        Dictionary containing the number of atoms of each 1-indexed fragment.
        For He--HOOH--Me cluster, `{1: 1, 2: 4, 3: 5}`.
    fragment_slice_dict
        Dictionary containing slices that index the gradient matrix for each of the 1-indexed fragments.
        For He--HOOH--Me cluster, `{1: slice(0, 1), 2: slice(1, 5), 3: slice(5, 10)}`.
    reverse
        If True, contract *hess* from supersystem size and shape that which is natural for *bas*.

    Returns
    -------
        Hessian array padded with zeros to supersystem size, (3 * <nat of supersystem>, 3 * <nat of supersystem>).
        If reverse=True, Hessian array extracted to natural size, (3 * <nat in bas>, 3 * <nat in bas>).

    """
    if reverse:
        nat = sum(fragment_size_dict[ifr] for ifr in bas)
    else:
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
            if reverse:
                ret[rel_sl1, rel_sl2] = hess[abs_sl1, abs_sl2]
            else:
                ret[abs_sl1, abs_sl2] = hess[rel_sl1, rel_sl2]

    return ret


def labeler(mc_level_lbl: str, frag: Tuple[int, ...], bas: Tuple[int, ...]) -> str:
    """Form label from model chemistry id and fragment and basis indices.

    Parameters
    ----------
    mc_level_lbl
        Key identifying the model chemistry. May be `"(auto)"`. Often the
        ManyBodyInput.specification.specification keys.
    frag
        List of 1-indexed fragments active in the supersystem.
    bas
        List of 1-indexed fragments with active basis sets in the supersystem.
        All those in *frag* plus any ghost.

    Returns
    -------
    str
        JSON string from inputs: `labeler("mp2",(1), (1, 2))` returns `'["mp2", 1, [1, 2]]'`.
    """
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


def collect_vars(
        bsse: str,
        prop: str,
        body_dict: Mapping[int, Union[float, np.ndarray]],
        max_nbody: int,
        embedding: bool = False,
        supersystem_ie_only: bool = False,
    ) -> Dict:
    """From *body_dict*, construct QCVariables.

    Parameters
    ----------
    bsse
        Uppercase label for a single many-body treatment, generally a value of BsseEnum.
    prop
        Uppercase label for a single property, generally a value of DriverEnum.
    body_dict
        Dictionary of minimal per-body info already specialized for property *prop* and treatment *bsse*. May contain either total data or interaction data, both cummulative not additive, from 1-body to max_nbody-body (see also *supersystem_ie_only*). Interaction data signaled by zero float or array for 1-body. May contain multiple model chemistry levels.
    max_nbody
        _description_
    embedding, optional
        Is charge embedding enabled, by default False?
    supersystem_ie_only, optional
        Is data available in *body_dict* only for 1-body (possibly zero) and nfr-body levels? By default False: data is available for consecutive levels, up to max_nbody-body.

    Returns
    -------
        _description_. Empty return if *embedding* enabled.
    """
    previous_e = body_dict[1]
    property_shape = find_shape(previous_e)
    tot_e = bool(np.count_nonzero(previous_e))
    nbody_range = list(body_dict)
    nbody_range.sort()
    res = {}
    if embedding:
        return res

    if tot_e:
        res[f"{bsse}-CORRECTED TOTAL {prop}"] = body_dict[max_nbody]
    res[f"{bsse}-CORRECTED INTERACTION {prop}"] = body_dict[max_nbody] - body_dict[1]
    res[f"{bsse}-CORRECTED INTERACTION {prop} THROUGH 1-BODY"] = shaped_zero(property_shape)

    if supersystem_ie_only:
        nfr = nbody_range[-1]
        for nb in [nfr]:
            res[f"{bsse}-CORRECTED INTERACTION {prop} THROUGH {nb}-BODY"] = body_dict[nb] - body_dict[1]
            if nb == 2:
                res[f"{bsse}-CORRECTED {nb}-BODY CONTRIBUTION TO {prop}"] = body_dict[nb] - body_dict[nb - 1]
        if tot_e:
            for nb in [1, nfr]:
                res[f"{bsse}-CORRECTED TOTAL {prop} THROUGH {nb}-BODY"] = body_dict[nb]
    else:
        for nb in range(2, max(nbody_range) + 1):
            res[f"{bsse}-CORRECTED INTERACTION {prop} THROUGH {nb}-BODY"] = body_dict[nb] - body_dict[1]
            res[f"{bsse}-CORRECTED {nb}-BODY CONTRIBUTION TO {prop}"] = body_dict[nb] - body_dict[nb - 1]
        if tot_e:
            for nb in nbody_range:
                res[f"{bsse}-CORRECTED TOTAL {prop} THROUGH {nb}-BODY"] = body_dict[nb]

    return res


def provenance_stamp(routine: str) -> Dict[str, str]:
    """Return dictionary satisfying QCSchema,
    https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/definitions.py#L23-L41
    with QCManyBody's credentials for creator and version. The
    generating routine's name is passed in through `routine`.

    """
    import qcmanybody
    return {"creator": "QCManyBody", "version": qcmanybody.__version__, "routine": routine}
