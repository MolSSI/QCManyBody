from importlib.metadata import version

from .core import ManyBodyCore
from .core import ManyBodyCalculator  # legacy near-alias to ManyBodyCore
from .computer import ManyBodyComputer
from .models import BsseEnum
from .utils import labeler, delabeler, resize_gradient, resize_hessian

__version__ = version("qcmanybody")
