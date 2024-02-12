from __future__ import annotations

from enum import Enum
from typing import Tuple


class BsseEnum(str, Enum):
    """Available basis-set superposition error (BSSE) treatments."""

    nocp = "nocp"  # plain supramolecular interaction energy
    cp = "cp"  # counterpoise correction
    vmfc = "vmfc"  # Valiron-Mayer function counterpoise


FragBasIndex = Tuple[Tuple[int], Tuple[int]]
