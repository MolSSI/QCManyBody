from importlib.metadata import version

__version__ = version("qcmanybody")
del version

# isort: off
from .core import ManyBodyCalculator  # legacy near-alias to ManyBodyCore
from .core import ManyBodyCore
from .computer import ManyBodyComputer
from .parallel import ParallelManyBodyExecutor, ParallelConfig

# isort: on
from .models import BsseEnum
from .utils import delabeler, labeler, resize_gradient, resize_hessian
