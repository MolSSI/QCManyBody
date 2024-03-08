from __future__ import annotations

from typing import Dict, Union, Set

import numpy as np

from .models import FragBasIndex
from .utils import labeler, zeros_like, all_same_shape


def sum_cluster_data(
    data: Dict[str, Union[float, np.ndarray]],
    compute_list: Set[FragBasIndex],
    mc_level_lbl: str,
    vmfc: bool = False,
    nb: int = 0,  # used only for vmfc
) -> Union[float, np.ndarray]:
    sign = 1

    # Check the shapes of the data. Should all be the same
    if not all_same_shape(data.values()):
        raise ValueError("All values in data dictionary must have the same shape.")

    first_key = next(iter(data))
    ret = zeros_like(data[first_key])

    for frag, bas in compute_list:
        ene = data[labeler(mc_level_lbl, frag, bas)]

        if vmfc:
            sign = (-1) ** (nb - len(frag))

        ret += sign * ene

    return ret
