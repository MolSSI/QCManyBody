from importlib.metadata import version

from .manybody import ManyBodyCalculator
from .models import BsseEnum
from .utils import labeler, delabeler

__version__ = version("qcmanybody")
