import pytest
import qcelemental
from qcelemental import constants
from qcelemental.models import Molecule
from qcelemental.testing import compare_recursive

from qcmanybody import ManyBodyCalculator  # test old name still operational
from qcmanybody import ManyBodyComputer
from qcmanybody.models import AtomicSpecification, ManyBodyInput


@pytest.fixture
def he_tetramer():
    a2 = 2 / constants.bohr2angstroms
    return Molecule(symbols=["He", "He", "He", "He"], fragments=[[0], [1], [2], [3]], geometry=[0, 0, 0, 0, 0, a2, 0, a2, 0, 0, a2, a2])


def test_noncontiguous_fragments_ordinary():
    with pytest.raises(qcelemental.exceptions.ValidationError) as e:
        Molecule(symbols=["H", "Ne", "Cl"], geometry=[0, 0, 0, 2, 0, 0, 0, 2, 0], fragments=[[0, 2], [1]])

    assert "QCElemental would need to reorder atoms to accommodate non-contiguous fragments" in str(e.value)


def test_noncontiguous_fragments_evaded():
    neon_in_hcl = Molecule(symbols=["H", "Ne", "Cl"], geometry=[0, 0, 0, 2, 0, 0, 0, 2, 0], fragments=[[0, 2], [1]], validated=True)

    with pytest.raises(ValueError) as e:
        ManyBodyCalculator(neon_in_hcl, ["cp"], {2: "mp2", 1: "mp2"}, False, False, None)

    assert "QCManyBody: non-contiguous fragments could be implemented but aren't at present" in str(e.value)


def test_noncontiguous_nbody_levels_same_mc(he_tetramer):
    with pytest.raises(ValueError) as e:
        ManyBodyCalculator(he_tetramer, ["cp"], {2: "mp2", 1: "mp2", 4: "mp2", 3: "ccsd"}, True, False, None)

    assert "QCManyBody: N-Body levels must be contiguous within a model chemistry spec" in str(e.value)


@pytest.mark.parametrize("mbe_keywords,ref_count,ref_text", [
    pytest.param(
        {"bsse_type": ["nocp", "cp", "vmfc"]},
        # 65,
        {"all": {"(auto)": {4: 1, 3: 8, 2: 24, 1: 32}},
         "cp": {"(auto)": {4: 1, 3: 4, 2: 6, 1: 4}},  # other 4 1b req'd are in nocp
         "nocp": {"(auto)": {4: 1, 3: 4, 2: 6, 1: 4}},
         "vmfc_compute": {"(auto)": {4: 1, 3: 8, 2: 24, 1: 32}},
         },
        """
    Model chemistry "(auto)" (§A):         65
        Number of 1-body computations:     32 (nocp: 4, cp: 4, vmfc_compute: 32)
        Number of 2-body computations:     24 (nocp: 6, cp: 6, vmfc_compute: 24)
        Number of 3-body computations:      8 (nocp: 4, cp: 4, vmfc_compute: 8)
        Number of 4-body computations:      1 (nocp: 1, cp: 1, vmfc_compute: 1)""",
        id="4b_all"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 3},
        # 18,
        {"all": {"(auto)": {1: 8, 2: 6, 3: 4}},
         "cp": {"(auto)": {1: 4, 2: 6, 3: 4}},
         "nocp": {"(auto)": {1: 4}},
         "vmfc_compute": {"(auto)": {}},
        },
        """
    Model chemistry "(auto)" (§A):         18
        Number of 1-body computations:      8 (nocp: 4, cp: 4, vmfc_compute: 0)
        Number of 2-body computations:      6 (nocp: 0, cp: 6, vmfc_compute: 0)
        Number of 3-body computations:      4 (nocp: 0, cp: 4, vmfc_compute: 0)""",
        id="3b_cp_rtd"),
    pytest.param(
        {"bsse_type": "vmfc", "return_total_data": False, "max_nbody": 3},
        # 50,
        {"all": {"(auto)": {1: 28, 2: 18, 3: 4}},
         "cp": {"(auto)": {}},
         "nocp": {"(auto)": {1: 4, 2: 6, 3: 4}},  # free with vmfc
         "vmfc_compute": {"(auto)": {1: 28, 2: 18, 3: 4}},
        },
        """
    Model chemistry "(auto)" (§A):         50
        Number of 1-body computations:     28 (nocp: 4, cp: 0, vmfc_compute: 28)
        Number of 2-body computations:     18 (nocp: 6, cp: 0, vmfc_compute: 18)
        Number of 3-body computations:      4 (nocp: 4, cp: 0, vmfc_compute: 4)""",
        id="3b_vmfc"),
])
def test_count_he4_single(mbe_keywords, ref_count, ref_text, he_tetramer):
    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": "mybas"}, program="myqc", driver="energy")
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    ret = ManyBodyComputer.from_manybodyinput(mbe_model, build_tasks=False)
    ret = ret.qcmb_core

    text, dcount = ret.format_calc_plan()
    print(text)
    assert compare_recursive(ref_count["all"], dcount, atol=1.e-6)
    for sset in ["all", "cp", "nocp", "vmfc_compute"]:
        text, dcount = ret.format_calc_plan(sset)
        assert compare_recursive(ref_count[sset], dcount, atol=1.e-6)
    assert ref_text.strip() == text.strip(), f"Output printing failed\n{text}\n{ref_text}"


@pytest.mark.parametrize("mbe_keywords,ref_count,ref_text", [
    pytest.param(
        {"bsse_type": ["nocp", "cp", "vmfc"], "levels": {4: "c4-mp2", 2: "c4-ccsd"},},
        {"all": {"c4-ccsd": {2: 12, 1: 20}, "c4-mp2": {4: 1, 3: 8, 2: 24, 1: 20}},
         "cp": {"c4-ccsd": {2: 6, 1: 4}, "c4-mp2": {4: 1, 3: 4, 2: 6, 1: 4}},
         "nocp": {"c4-ccsd": {2: 6, 1: 4}, "c4-mp2": {4: 1, 3: 4, 2: 6, 1: 4}},
         "vmfc_compute": {"c4-ccsd": {2: 6, 1: 16}, "c4-mp2": {4: 1, 3: 8, 2: 18, 1: 16}},
         },
         """
    Model chemistry "c4-ccsd" (§A):        32
        Number of 1-body computations:     20 (nocp: 4, cp: 4, vmfc_compute: 16)
        Number of 2-body computations:     12 (nocp: 6, cp: 6, vmfc_compute: 6)

    Model chemistry "c4-mp2" (§B):         53
        Number of 1-body computations:     20 (nocp: 4, cp: 4, vmfc_compute: 16)
        Number of 2-body computations:     24 (nocp: 6, cp: 6, vmfc_compute: 18)
        Number of 3-body computations:      8 (nocp: 4, cp: 4, vmfc_compute: 8)
        Number of 4-body computations:      1 (nocp: 1, cp: 1, vmfc_compute: 1)""",
        id="4b_all"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 3, "levels": {3: "c4-mp2", 2: "c4-ccsd"},},
        {"all": {"c4-ccsd": {1: 8, 2: 6}, "c4-mp2": {1: 4, 2: 6, 3: 4}},
         "cp": {"c4-ccsd": {1: 4, 2: 6}, "c4-mp2": {1: 4, 2: 6, 3: 4}},
         "nocp": {"c4-ccsd": {1: 4}, "c4-mp2": {}},
         "vmfc_compute": {"c4-ccsd": {}, "c4-mp2": {}},
        },
        """
    Model chemistry "c4-ccsd" (§A):        14
        Number of 1-body computations:      8 (nocp: 4, cp: 4, vmfc_compute: 0)
        Number of 2-body computations:      6 (nocp: 0, cp: 6, vmfc_compute: 0)

    Model chemistry "c4-mp2" (§B):         14
        Number of 1-body computations:      4 (nocp: 0, cp: 4, vmfc_compute: 0)
        Number of 2-body computations:      6 (nocp: 0, cp: 6, vmfc_compute: 0)
        Number of 3-body computations:      4 (nocp: 0, cp: 4, vmfc_compute: 0)""",
        id="3b_cp_rtd"),
   pytest.param(
       {"bsse_type": "vmfc", "return_total_data": False, "max_nbody": 3, "levels": {2: "c4-ccsd", 3: "ac4-mp2"},},
       {"all": {"c4-ccsd": {2: 6, 1: 16}, "ac4-mp2": {3: 4, 2: 12, 1: 12}},
         "cp": {"c4-ccsd": {}, "ac4-mp2": {}},
         "nocp": {"c4-ccsd": {}, "ac4-mp2": {}},  # no free nocp with vmfc for multilevel
         "vmfc_compute": {"c4-ccsd": {2: 6, 1: 16}, "ac4-mp2": {3: 4, 2: 12, 1: 12}},
       },
       """
    Model chemistry "c4-ccsd" (§A):        22
        Number of 1-body computations:     16 (nocp: 0, cp: 0, vmfc_compute: 16)
        Number of 2-body computations:      6 (nocp: 0, cp: 0, vmfc_compute: 6)

    Model chemistry "ac4-mp2" (§B):        28
        Number of 1-body computations:     12 (nocp: 0, cp: 0, vmfc_compute: 12)
        Number of 2-body computations:     12 (nocp: 0, cp: 0, vmfc_compute: 12)
        Number of 3-body computations:      4 (nocp: 0, cp: 0, vmfc_compute: 4)""",
       id="3b_vmfc"),
    pytest.param(
        {"bsse_type": "nocp", "levels": {2: "c4-ccsd", "supersystem": "mp3"}},
        {"nocp": {"c4-ccsd": {1: 4, 2: 6}, "mp3": {2: 6, 1: 4, 4: 1}}, "cp": {"c4-ccsd": {}, "mp3": {}}, "vmfc_compute": {"c4-ccsd": {}, "mp3": {}}, "all": {"c4-ccsd": {1: 4, 2: 6}, "mp3": {2: 6, 1: 4, 4: 1}}},
        """
    Model chemistry "c4-ccsd" (§A):        10
        Number of 1-body computations:      4 (nocp: 4, cp: 0, vmfc_compute: 0)
        Number of 2-body computations:      6 (nocp: 6, cp: 0, vmfc_compute: 0)

    Model chemistry "mp3" (§B):            11
        Number of 1-body computations:      4 (nocp: 4, cp: 0, vmfc_compute: 0)
        Number of 2-body computations:      6 (nocp: 6, cp: 0, vmfc_compute: 0)
        Number of 4-body computations:      1 (nocp: 1, cp: 0, vmfc_compute: 0)""",
        id="2b_uncp_ss"
    ),
])
def test_count_he4_multi(mbe_keywords, ref_count, ref_text, he_tetramer, request):
    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": "mybas"}, program="myqc", driver="energy")
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    ret = ManyBodyComputer.from_manybodyinput(mbe_model, build_tasks=False)
    ret = ret.qcmb_core

    text, dcount = ret.format_calc_plan()
    print(text)
    assert compare_recursive(ref_count["all"], dcount, atol=1.e-6)
    for sset in ["all", "cp", "nocp", "vmfc_compute"]:
        text, dcount = ret.format_calc_plan(sset)
        assert compare_recursive(ref_count[sset], dcount, atol=1.e-6)
    assert ref_text.strip() == text.strip(), f"Output printing failed\n{text}\n{ref_text}"
