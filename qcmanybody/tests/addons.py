from typing import List

import pytest

from qcelemental.util import parse_version, which, which_import
from qcengine.testing import _programs as _programs_qcng


__all__ = [
    "using",
    "uusing",
]


def is_qcfractal_new_enough(version_feature_introduced):
    if not which_import('qcfractal', return_bool=True):
        return False
    import qcfractal
    return parse_version(qcfractal.__version__) >= parse_version(version_feature_introduced)


# Figure out what is imported
# * anything _not_ in QCEngine goes here
_programs = {
    "vasp": False,
}


def has_program(name):
    # any aliases or merged names go here
    if name in _programs:
        return _programs[name]
    elif name in _programs_qcng:
        return _programs_qcng[name]
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
