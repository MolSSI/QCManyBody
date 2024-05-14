from importlib.metadata import version

from .manybody import ManyBodyCalculator
from .models import BsseEnum
from .utils import labeler, delabeler, resize_gradient, resize_hessian
from .qcng_computer import ManyBodyComputerQCNG

__version__ = version("qcmanybody")
