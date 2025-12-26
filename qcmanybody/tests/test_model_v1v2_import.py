import sys
import warnings

import pytest


def _test_instantiate_mdl(clsins):
    mbe_model = {"specification": {"specification": {"model": {"method": "mp2", "basis": "cc-pvdz"}, "driver": "energy", "program": "myqc"}, "driver": "energy", "keywords": {"bsse_type": "cp"}}, "molecule": {"geometry": [0, 0, 0, 2, 0, 0], "symbols": ["he", "he"]}}

    if sys.version_info < (3, 14):
        clsins(**mbe_model)

    else:
        with pytest.raises(RuntimeError) as exc:
            clsins(**mbe_model)

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
