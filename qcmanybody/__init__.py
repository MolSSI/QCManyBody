from importlib.metadata import version

# isort: off
from .core import ManyBodyCalculator  # legacy near-alias to ManyBodyCore
from .core import ManyBodyCore
from .computer import ManyBodyComputer

# isort: on
from .models import BsseEnum
from .utils import delabeler, labeler, resize_gradient, resize_hessian

__version__ = version("qcmanybody")
del version
