"""
Safety shim for QCSchema v1 imports.

This module intentionally avoids importing v1 model internals eagerly. If the
Pydantic v1 compatibility layer (``pydantic.v1``) is useable, we will try to
re-export the real classes. Otherwise this module exposes placeholder classes
that emit a deprecation warning on import and raise a clear RuntimeError when
an attempt is made to instantiate any v1 model.
"""

from __future__ import annotations

import importlib
import sys
import warnings

_MSG = (
    "qcmanybody.models.v1 is active but incompatible with Python 3.14+ "
    "(pydantic.v1 is not available). Imports will provide non-functional placeholders; "
    "instantiating any v1 model will raise a RuntimeError. Please run on Python < 3.14 "
    "or migrate to qcmanybody.models.v2."
)

# Warn on import so users see the incompatibility when they explicitly import v1
# Use FutureWarning (visible by default) so users notice the issue during import
warnings.warn(_MSG, FutureWarning, stacklevel=2)


def _make_placeholder(name: str):
    """Create a placeholder class for a v1 model.

    The class is importable (so `from ... import Name` works) but raises a
    RuntimeError on instantiation with an actionable message.
    """
    _MSG2 = (
        f"QCSchema v1 model '{name}' cannot be instantiated in this environment. "
        + "Reason: pydantic.v1 is unavailable on Python 3.14+. "
        + "Use qcmanybody.models.v2 or run Python <3.14. See docs/MIGRATION.md"
    )

    def __init__(self, *args, **kwargs):
        raise RuntimeError(_MSG2)

    def __repr__(self):
        return f"<Unavailable QCSchema v1 model {name}>"

    return type(name, (), {"__init__": __init__, "__repr__": __repr__})


# Names this module should export (keeps parity with the previous file layout)
_EXPORT_NAMES = [
    "AtomicSpecification",
    "BsseEnum",
    "FragBasIndex",
    "ManyBodyInput",
    "ManyBodyKeywords",
    "ManyBodyProtocols",
    "ManyBodySpecification",
    "MAX_NBODY",
    "ManyBodyResult",
    "ManyBodyResultProperties",
]


def _use_real_if_possible():
    """Attempt to import and re-export real v1 symbols when pydantic.v1 is present.

    Falls back silently to placeholders on any error.
    """
    # Note: can't test on `import pydantic.v1` b/c it's not necessarily functional
    if sys.version_info >= (3, 14):
        return False

    # Map where names are defined in the original layout. Import cautiously.
    mapping = {
        "AtomicSpecification": (".manybody_input_pydv1", "AtomicSpecification"),
        "BsseEnum": (".manybody_input_pydv1", "BsseEnum"),
        "FragBasIndex": (".manybody_input_pydv1", "FragBasIndex"),
        "ManyBodyInput": (".manybody_input_pydv1", "ManyBodyInput"),
        "ManyBodyKeywords": (".manybody_input_pydv1", "ManyBodyKeywords"),
        "ManyBodyProtocols": (".manybody_input_pydv1", "ManyBodyProtocols"),
        "ManyBodySpecification": (".manybody_input_pydv1", "ManyBodySpecification"),
        "MAX_NBODY": (".manybody_output_pydv1", "MAX_NBODY"),
        "ManyBodyResults": (".manybody_output_pydv1", "ManyBodyResult"),
        "ManyBodyResultProperties": (".manybody_output_pydv1", "ManyBodyResultProperties"),
    }

    pkg = __name__.rsplit(".", 1)[0]

    for name, (submod, attr) in mapping.items():
        try:
            module = importlib.import_module(submod, package="qcmanybody.models.v1")
            value = getattr(module, attr)
            globals()[name] = value
        except Exception:
            # Leave placeholder if anything goes wrong
            globals()[name] = _make_placeholder(name)

    return True


# First create placeholders for all names so imports always succeed
for _n in _EXPORT_NAMES:
    globals()[_n] = _make_placeholder(_n)


# Then attempt to shadow them with real implementations when available
_use_real_if_possible()
