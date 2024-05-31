from importlib.metadata import version

from .manybody import ManyBodyCore
from .manybody import ManyBodyCore as ManyBodyCalculator  # legacy alias
from .computer import ManyBodyComputer
from .models import BsseEnum
from .utils import labeler, delabeler, resize_gradient, resize_hessian

__version__ = version("qcmanybody")
