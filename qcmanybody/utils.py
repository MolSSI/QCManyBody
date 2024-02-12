from __future__ import annotations

import json
import numpy as np
from typing import Tuple, Dict
from qcelemental import constants


def labeler(mc_level_lbl: str, frag: Tuple[int, ...], bas: Tuple[int, ...]) -> str:
    return str(mc_level_lbl) + "_" + json.dumps((frag, bas))


def delabeler(item: str) -> Tuple[str, Tuple[int, ...], Tuple[int, ...]]:
    """Transform labels like string "1_((2,), (1, 2))" into tuple (1, (2,), (1, 2))."""

    mc, _, fragbas = item.partition("_")
    frag, bas = json.loads(fragbas)
    return str(mc), frag, bas


def shaped_zero(der: int, nat: int):
    if der == 0:
        return 0.0
    elif der == 1:
        arr_shape = (nat, 3)
        return np.zeros(arr_shape)
    elif der == 2:
        arr_shape = (nat * 3, nat * 3)
        return np.zeros(arr_shape)
    else:
        raise ValueError(f"Invalid derivative level {der}")


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
