from typing import List

import pytest
from qcelemental.util import parse_version, which, which_import

try:
    from qcengine.testing import _programs as _programs_qcng
except ModuleNotFoundError:
    _programs_qcng = {}


__all__ = [
    "using",
    "uusing",
]


def is_qcfractal_new_enough(version_feature_introduced):
    if not which_import("qcfractal", return_bool=True):
        return False
    import qcfractal

    return parse_version(qcfractal.__version__) >= parse_version(version_feature_introduced)


# Figure out what is imported
# * anything _not_ in QCEngine goes here
# * to allow qcng as optional dependency, duplicate the using's in qcmb's test suite here
_programs = {
    "vasp": False,
    "cfour": which("xcfour", return_bool=True),
    "gamess": which("rungms", return_bool=True),
    "geometric": which_import("geometric", return_bool=True),
    "optking": which_import("optking", return_bool=True),
    "nwchem": which("nwchem", return_bool=True),
    "psi4": which("psi4", return_bool=True),
    "qcengine": which_import("qcengine", return_bool=True),
}


def has_program(name):
    # any aliases or merged names go here
    if name in _programs:
        return _programs[name]
    elif name in _programs_qcng:
        return _programs_qcng[name]
    elif name in ["optking_genopt", "geometric_genopt"]:
        # give up rather than duplicate
        return False
    else:
        raise KeyError(f"Program {name} not registered with QCManyBody testing.")


_using_cache = {}


def _using(program: str) -> None:
    if program not in _using_cache:
        import_message = f"Not detecting module {program}. Install package if necessary to enable tests."
        skip = pytest.mark.skipif(has_program(program) is False, reason=import_message)
        general = pytest.mark.addon
        particular = getattr(pytest.mark, program)

        all_marks = (skip, general, particular)
        _using_cache[program] = [_compose_decos(all_marks), all_marks]


def _compose_decos(decos):
    # thanks, https://stackoverflow.com/a/45517876
    def composition(func):
        for deco in reversed(decos):
            func = deco(func)
        return func

    return composition


def uusing(program: str):
    """Apply 3 marks: skipif program not detected, label "addon", and label program.
    This is the decorator form for whole test functions: `@mark\n@mark`.

    """
    _using(program)
    return _using_cache[program][0]


def using(program: str) -> List:
    """Apply 3 marks: skipif program not detected, label "addon", and label program.
    This is the inline form for parameterizations: `marks=[]`.
    In combo, do `marks=[*using(), pytest.mark.quick]`

    """
    _using(program)
    return _using_cache[program][1]


@pytest.fixture(scope="function", params=[None, "v1", "v2"])
def schema_versions(request):
    import qcmanybody

    if request.param == "v1":
        return qcmanybody.models.v1
    elif request.param == "v2":
        try:
            import qcmanybody.models.v2
        except (ModuleNotFoundError, ImportError):
            pytest.skip("QCManyBody v2 not tested since QCElemental v2 not available.")
        return qcmanybody.models.v2
    else:
        return qcmanybody.models
