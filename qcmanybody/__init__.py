from importlib.metadata import version

__version__ = version("qcmanybody")
del version

# isort: off
from .core import ManyBodyCalculator  # legacy near-alias to ManyBodyCore
from .core import ManyBodyCore
from .computer import ManyBodyComputer

# isort: on
# BsseEnum used mostly for values, not object, so ensure available. v1 to minimize surprise.
import sys

if sys.version_info < (3, 14):
    from .models.v1 import BsseEnum
else:
    from .models.v2 import BsseEnum
del sys

from .utils import delabeler, labeler, resize_gradient, resize_hessian
