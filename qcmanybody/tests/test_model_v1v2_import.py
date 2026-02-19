import sys
import warnings

import pytest

_mbe_model = {"molecule": {"geometry": [0, 0, 0, 2, 0, 0], "symbols": ["he", "he"], "fragments": [[0], [1]]}, "specification": {"driver": "energy", "keywords": {"bsse_type": "cp"}, "specification": {"model": {"method": "mp2", "basis": "cc-pvdz"}, "driver": "energy", "program": "myqc"}}}

def _test_instantiate_mdl(clsins):
    if sys.version_info < (3, 14):
        clsins(**_mbe_model)

    else:
        with pytest.raises(RuntimeError) as exc:
            clsins(**_mbe_model)

        assert "pydantic.v1 is unavailable" in str(exc.value)


def test_ins_mdl_var_models_models():
    from qcmanybody.models import ManyBodyInput as MyMBIn

    _test_instantiate_mdl(MyMBIn)


def test_ins_mdl_var_models_models_v1():
    from qcmanybody.models.v1 import ManyBodyInput as MyMBIn

    _test_instantiate_mdl(MyMBIn)


def test_ins_mdl_var_models_models_input_pydv1():
    from qcmanybody.models.manybody_input_pydv1 import ManyBodyInput as MyMBIn

    _test_instantiate_mdl(MyMBIn)


def _test_instantiate_cptr(clsins):

    if sys.version_info < (3, 14):
        clsins.from_manybodyinput(_mbe_model, build_tasks=False)

    else:
        with pytest.raises(RuntimeError) as exc:
            clsins.from_manybodyinput(_mbe_model, build_tasks=False)

        assert "pydantic.v1 is unavailable" in str(exc.value)


def test_ins_mdl_var_computer():
    from qcmanybody import ManyBodyComputer as MyComp

    _test_instantiate_cptr(MyComp)


def test_ins_mdl_var_computer_v1():
    from qcmanybody.v1 import ManyBodyComputer as MyComp

    _test_instantiate_cptr(MyComp)


def test_ins_mdl_var_computer_v1_computer_error():

    if sys.version_info < (3, 14):
        from qcmanybody.v1.computer import ManyBodyComputer as MyComp

        MyComp.from_manybodyinput(_mbe_model, build_tasks=False)

    else:
        # not a recc. way to import
        #   this can't get beyond the description=...__fields__ part of instantiation
        with pytest.raises(AttributeError) as exc:
            from qcmanybody.v1.computer import ManyBodyComputer as MyComp

        assert "has no attribute '__fields__'" in str(exc.value)
