import sys

try:
    import pydantic
except ImportError:  # pragma: no cover
    raise ImportError(
        "Python module pydantic not found. Solve by installing it: "
        "`conda install pydantic -c conda-forge` or `pip install pydantic`"
    )

from . import v1

# someday, `from . import v1, v2`. for now, this lets qcmb work with pre-next (v0.50) QCElemental
try:
    from . import v2
except (ModuleNotFoundError, ImportError):
    pass

from .v1 import *

# Note that changing .v1 as default requires changing the shim classes in this dir, too.
