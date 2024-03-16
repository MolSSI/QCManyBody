"""
Tests the DQM compute dispatch module
"""

import pytest
import pydantic

from qcelemental.models import DriverEnum, Molecule
from qcmanybody.models.manybody_v1 import BsseEnum, ManyBodyInput

import qcengine as qcng

# qcng: from qcengine.procedures.manybody import ManyBodyComputerQCNG
from qcmanybody.models.qcng_computer import ManyBodyComputerQCNG


@pytest.fixture(scope="function")
def mbe_data():
    henehh = Molecule(symbols=["He", "Ne", "H", "H"], fragments=[[0], [1], [2, 3]], geometry=[0, 0, 0, 2, 0, 0, 0, 1, 0, 0, -1, 0])
    return {
        "specification": {
            "specification": {
                "model": {
                    "method": "hf",
                    "basis": "cc-pvdz",
                },
                "driver": "energy",
                "program": "psi4",
            },
            "keywords": {
                "bsse_type": "cp",
            },
            "driver": "gradient",
        },
        "molecule": henehh,
    }


@pytest.fixture(scope="function")
def mbe_data_multilevel():
    henehh = Molecule(symbols=["He", "Ne", "H", "H"], fragments=[[0], [1], [2, 3]], geometry=[0, 0, 0, 2, 0, 0, 0, 1, 0, 0, -1, 0])
    return {
        "specification": {
            "specification": {
                "(auto)": { #p4-hf-dz": {
                    "model": {
                        "method": "hf",
                        "basis": "cc-pvdz",
                    },
                    "driver": "energy",
                    "program": "psi4",
                },
                "p4-mp2-dz": {
                    "model": {
                        "method": "mp2",
                        "basis": "cc-pvdz",
                    },
                    "driver": "energy",
                    "program": "psi4",
                },
                "c4-hf-dz": {
                    "model": {
                        "method": "hf",
                        "basis": "cc-pvdz",
                    },
                    "driver": "energy",
                    "program": "cfour",
                },
                "c4-ccsd-tz": {
                    "model": {
                        "method": "ccsd",
                        "basis": "cc-pvtz",
                    },
                    "driver": "energy",
                    "program": "cfour",
                },
            },
            "keywords": {
                "bsse_type": "nocp",
            },
            "driver": "gradient",
        },
        "molecule": henehh,
    }


@pytest.mark.parametrize("driver,kws,ans", [
    pytest.param("energy", {}, False),
    pytest.param("energy", {"return_total_data": True}, True),
    pytest.param("energy", {"return_total_data": False}, False),
    pytest.param("gradient", {}, True),
    pytest.param("hessian", {}, True),
    pytest.param("hessian", {"return_total_data": True}, True),
    pytest.param("hessian", {"return_total_data": False}, False),
])
def test_mbe_rtd(mbe_data, driver, kws, ans):
    mbe_data["specification"]["driver"] = driver
    mbe_data["specification"]["keywords"] = kws

    input_model = ManyBodyInput(**mbe_data)
    comp_model = ManyBodyComputerQCNG.from_qcschema(input_model)

    assert comp_model.driver == driver
    assert comp_model.return_total_data == ans


@pytest.mark.parametrize("kws,ans", [
    pytest.param({}, [3, {3: "(auto)"}, [[1, 2, 3]] ]),
    pytest.param({"max_nbody": 3}, [3, {3: "(auto)"}, [[1, 2, 3]] ]),
    pytest.param({"max_nbody": 2}, [2, {2: "(auto)"}, [[1, 2]] ]),
    pytest.param({"max_nbody": 1}, [1, {1: "(auto)"}, [[1]] ]),

    # TODO? when supersystem_ie_only=T, nbodies_per_mc_level and levels isn't really accurate
    pytest.param({"supersystem_ie_only": True}, [3, {3: "(auto)"}, [[1, 2, 3]] ]),
    pytest.param({"supersystem_ie_only": False, "max_nbody": 2}, [2, {2: "(auto)"}, [[1, 2]] ]),
    pytest.param({"supersystem_ie_only": True, "max_nbody": 3}, [3, {3: "(auto)"}, [[1, 2, 3]] ]),

    pytest.param({"levels": {3: "mp2"}}, [3, {3: "mp2"}, [[1, 2, 3]] ]),
    pytest.param({"levels": {3: "ccsd", 2: "ccsd"}}, [3, {2: "ccsd", 3: "ccsd"}, [[1, 2], [3]] ]),
    pytest.param({"levels": {1: "mp2", 3: "ccsd"}}, [3, {1: "mp2", 3: "ccsd"}, [[1], [2, 3]] ]),
    pytest.param({"levels": {2: "ccsd", 3: "mp2"}}, [3, {2: "ccsd", 3: "mp2"}, [[1, 2], [3]] ]),
    pytest.param({"levels": {2: "ccsd"}}, [2, {2: "ccsd"}, [[1, 2]] ]),
    pytest.param({"levels": {2: "ccsd", 1: "ccsd(t)"}}, [2, {1: "ccsd(t)", 2: "ccsd"}, [[1], [2]] ]),
    pytest.param({"levels": {1: "ccsd(t)"}}, [1, {1: "ccsd(t)"}, [[1]] ]),
    pytest.param({"levels": {"supersystem": "hf", 1: "ccsd(t)"}}, [1, {1: "ccsd(t)", "supersystem": "hf"}, [[1], ["supersystem"]] ]),
    pytest.param({"levels": {2: "ccsd(t)", "supersystem": "hf"}}, [2, {2: "ccsd(t)", "supersystem": "hf"}, [[1, 2], ["supersystem"]] ]),

    pytest.param({"max_nbody": 3, "levels": {3: "mp2"}}, [3, {3: "mp2"}, [[1, 2, 3]] ]),
    pytest.param({"max_nbody": 3, "levels": {3: "ccsd", 2: "ccsd"}}, [3, {2: "ccsd", 3: "ccsd"}, [[1, 2], [3]] ]),
    pytest.param({"max_nbody": 3, "levels": {1: "mp2", 3: "ccsd"}}, [3, {1: "mp2", 3: "ccsd"}, [[1], [2, 3]] ]),
    pytest.param({"max_nbody": 3, "levels": {2: "ccsd", 3: "mp2"}}, [3, {2: "ccsd", 3: "mp2"}, [[1, 2], [3]] ]),
    pytest.param({"max_nbody": 2, "levels": {2: "ccsd"}}, [2, {2: "ccsd"}, [[1, 2]] ]),
    pytest.param({"max_nbody": 2, "levels": {2: "ccsd", 1: "ccsd(t)"}}, [2, {1: "ccsd(t)", 2: "ccsd"}, [[1], [2]] ]),
    pytest.param({"max_nbody": 1, "levels": {1: "ccsd(t)"}}, [1, {1: "ccsd(t)"}, [[1]] ]),
    pytest.param({"max_nbody": 1, "levels": {"supersystem": "hf", 1: "ccsd(t)"}}, [1, {1: "ccsd(t)", "supersystem": "hf"}, [[1], ["supersystem"]] ]),
    pytest.param({"max_nbody": 2, "levels": {2: "ccsd(t)", "supersystem": "hf"}}, [2, {2: "ccsd(t)", "supersystem": "hf"}, [[1, 2], ["supersystem"]] ]),

    pytest.param({"levels": {2: 'scf/sto-3g', 1: 'mp2/sto-3g'}}, [2, {1: 'mp2/sto-3g', 2: 'scf/sto-3g'}, [[1], [2]] ]),
    pytest.param({"levels": {1: 'mp2/sto-3g', 'supersystem': 'scf/sto-3g'}}, [1, {1: 'mp2/sto-3g', 'supersystem': 'scf/sto-3g'}, [[1], ["supersystem"]] ]),
    pytest.param({"levels": {'supersystem': 'scf/sto-3g', 1: 'mp2/sto-3g'}}, [1, {1: 'mp2/sto-3g', 'supersystem': 'scf/sto-3g'}, [[1], ["supersystem"]] ]),

    #pytest.param({}, [ , {}, [] ]),
])
def test_mbe_level_bodies(mbe_data, kws, ans):
    mbe_data["specification"]["keywords"] = kws

    input_model = ManyBodyInput(**mbe_data)
    comp_model = ManyBodyComputerQCNG.from_qcschema(input_model)

    assert comp_model.nfragments == 3
    assert comp_model.max_nbody == ans[0]
    assert list(comp_model.levels.items()) == list(ans[1].items())  # compare as OrderedDict
    assert comp_model.nbodies_per_mc_level == ans[2]


@pytest.mark.parametrize("kws,ans", [
    pytest.param({}, [3, {3: "(auto)"}, [[1, 2, 3]], {"1_((1,), (1,))": ("hf", "psi4")}]),
    pytest.param({"levels": {3: "p4-mp2-dz"}}, [3, {3: "p4-mp2-dz"}, [[1, 2, 3]], {"1_((3,), (3,))": ("mp2", "psi4")} ]),
    pytest.param({"levels": {1: "p4-mp2-dz", 3: "c4-ccsd-tz"}}, [3, {1: "p4-mp2-dz", 3: "c4-ccsd-tz"}, [[1], [2, 3]], {"1_((1,), (1,))": ("mp2", "psi4"), "2_((1, 2, 3), (1, 2, 3))": ("ccsd", "cfour")} ]),
])
def hide_test_mbe_multilevel(mbe_data_multilevel, kws, ans):
    mbe_data_multilevel["specification"]["keywords"] = kws

    input_model = ManyBodyInput(**mbe_data_multilevel)
    comp_model = ManyBodyComputerQCNG.from_qcschema(input_model, build_tasks=True)

    assert comp_model.nfragments == 3
    assert comp_model.max_nbody == ans[0]
    assert list(comp_model.levels.items()) == list(ans[1].items())  # compare as OrderedDict
    assert comp_model.nbodies_per_mc_level == ans[2]

    import pprint
    pprint.pprint(comp_model.model_dump(), width=200)

    for k, v in ans[3].items():
        assert comp_model.task_list[k].method == v[0]
        assert comp_model.task_list[k].program == v[1]


@pytest.mark.parametrize("kws,ans", [
    pytest.param({"levels": {5: "hi"}}, [5, {5: "hi"}, [[1, 2, 3, 4, 5]] ]),
    pytest.param({"levels": {5: "hi", 4: "md"}}, [5, {4: "md", 5: "hi"}, [[1, 2, 3, 4], [5]] ]),
    pytest.param({"levels": {5: "hi", 3: "md"}}, [5, {3: "md", 5: "hi"}, [[1, 2, 3], [4, 5]] ]),
    pytest.param({"levels": {5: "hi", 2: "md"}}, [5, {2: "md", 5: "hi"}, [[1, 2], [3, 4, 5]] ]),
    pytest.param({"levels": {5: "hi", 1: "md"}}, [5, {1: "md", 5: "hi"}, [[1], [2, 3, 4, 5]] ]),

    pytest.param({"levels": {4: "md", 3: "lo", 5: "hi"}}, [5, {3: "lo", 4: "md", 5: "hi"}, [[1, 2, 3], [4], [5]] ]),
    pytest.param({"levels": {4: "md", 1: "lo", 5: "hi"}}, [5, {1: "lo", 4: "md", 5: "hi"}, [[1], [2, 3, 4], [5]] ]),
    pytest.param({"levels": {4: "md", 2: "lo", 5: "hi"}}, [5, {2: "lo", 4: "md", 5: "hi"}, [[1, 2], [3, 4], [5]] ]),
    pytest.param({"levels": {3: "md", 1: "lo", 5: "hi"}}, [5, {1: "lo", 3: "md", 5: "hi"}, [[1], [2, 3], [4, 5]] ]),

    pytest.param({"levels": {3: "md", 1: "lo", 4: "hi"}}, [4, {1: "lo", 3: "md", 4: "hi"}, [[1], [2, 3], [4]] ]),
])
def test_mbe_level_5mer(mbe_data, kws, ans):
    he3ne2 = Molecule(symbols=["He", "He", "He", "Ne", "Ne"], fragments=[[0], [1], [2], [3], [4]], geometry=[0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 2, 0, 0, -2, 0])
    mbe_data["molecule"] = he3ne2
    mbe_data["specification"]["keywords"] = kws

    input_model = ManyBodyInput(**mbe_data)
    comp_model = ManyBodyComputerQCNG.from_qcschema(input_model)

    assert comp_model.nfragments == 5
    assert comp_model.max_nbody == ans[0]
    assert list(comp_model.levels.items()) == list(ans[1].items())  # compare as OrderedDict
    assert comp_model.nbodies_per_mc_level == ans[2]


@pytest.mark.parametrize("kws,errmsg", [
    pytest.param({"max_nbody": -2}, "should be between 1 and 3"),
    pytest.param({"max_nbody": -1}, "should be between 1 and 3"),
    pytest.param({"max_nbody": 4}, "should be between 1 and 3"),
    # fails on v1 b/c val 2 coerced to a str  pytest.param({"levels": {1: 2, 3: "mp2", 2: "ccsd"}}, "asdf"),  # `1: 2 is old syntax and doesn't pass typing
    pytest.param({"max_nbody": 1, "supersystem_ie_only": True}, "Cannot skip intermediate n-body jobs")
    # [3, supersystem]
])
def test_mbe_level_fails(mbe_data, kws, errmsg):
    mbe_data["specification"]["keywords"] = kws

    # v2: with pytest.raises(Exception):
    with pytest.raises((pydantic.ValidationError, pydantic.v1.ValidationError)) as e:
        input_model = ManyBodyInput(**mbe_data)
        ManyBodyComputerQCNG.from_qcschema(input_model)

    assert errmsg in str(e.value)


@pytest.mark.parametrize("kws,ans", [
    pytest.param({}, [BsseEnum.cp]),
    pytest.param({"bsse_type": "CP"}, [BsseEnum.cp]),
    pytest.param({"bsse_type": "nocp"}, [BsseEnum.nocp]),
    pytest.param({"bsse_type": ["vmfc"]}, [BsseEnum.vmfc]),
    pytest.param({"bsse_type": ["vmfc", "nocp"]}, [BsseEnum.vmfc, BsseEnum.nocp]),
    pytest.param({"bsse_type": ["ssFC", "nocp"]}, [BsseEnum.cp, BsseEnum.nocp]),
    pytest.param({"bsse_type": ["ssfc", "cp"]}, [BsseEnum.cp]),
    pytest.param({"bsse_type": ["ssfc", "vmfc", "nocp", "cp"]}, [BsseEnum.cp, BsseEnum.vmfc, BsseEnum.nocp]),
    pytest.param({"bsse_type": "mycp"}, "error"),
    pytest.param({"bsse_type": ["CP", "mycp"]}, "error"),
])
def test_mbe_bsse_type(mbe_data, kws, ans):
    mbe_data["specification"]["keywords"] = kws

    if ans == "error":
        with pytest.raises((pydantic.ValidationError, pydantic.v1.ValidationError)) as e:
            input_model = ManyBodyInput(**mbe_data)

        assert "not a valid enumeration member; permitted: 'nocp', 'cp', 'vmfc'" in str(e.value)
        return

    input_model = ManyBodyInput(**mbe_data)
    comp_model = ManyBodyComputerQCNG.from_qcschema(input_model)

    assert comp_model.bsse_type == ans


@pytest.mark.parametrize("kws,ans", [
    pytest.param({}, False),
    pytest.param({"max_nbody": 3, "supersystem_ie_only": True}, True),
    pytest.param({"supersystem_ie_only": False}, False),
    pytest.param({"max_nbody": 1, "supersystem_ie_only": True},
        "Cannot skip intermediate n-body jobs when max_nbody=1 != nfragments"),
])
def test_mbe_sie(mbe_data, kws, ans):
    mbe_data["specification"]["keywords"] = kws

    if isinstance(ans, str):
        input_model = ManyBodyInput(**mbe_data)
        with pytest.raises((pydantic.ValidationError, pydantic.v1.ValidationError)) as e:
            comp_model = ManyBodyComputerQCNG.from_qcschema(input_model)

        assert ans in str(e.value)
        return

    input_model = ManyBodyInput(**mbe_data)
    comp_model = ManyBodyComputerQCNG.from_qcschema(input_model)

    assert comp_model.supersystem_ie_only == ans

