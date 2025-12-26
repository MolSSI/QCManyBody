import importlib
import sys
import warnings

_MSG = (
    "qcmanybody.v1 is active but incompatible with Python 3.14+ "
    "(pydantic.v1 is not available). Imports will provide non-functional placeholders; "
    "instantiating any v1 model will raise a RuntimeError. Please run on Python < 3.14 "
    "or migrate to qcmanybody.v2."
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
    "ManyBodyComputer",
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
        "ManyBodyComputer": (".computer", "ManyBodyComputer"),
    }

    pkg = __name__.rsplit(".", 1)[0]

    for name, (submod, attr) in mapping.items():
        try:
            module = importlib.import_module(submod, package="qcmanybody.v1")
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
