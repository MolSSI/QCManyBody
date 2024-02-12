from qcmanybody.models import BsseEnum
from qcmanybody.builder import build_nbody_compute_list


def test_build_nbody_compute_list_1():
    # Single fragment, only asking for 1 body CP and NoCP
    r = build_nbody_compute_list(
        bsse_type=[BsseEnum.nocp, BsseEnum.cp],
        nfragments=1,
        nbodies=[1],
        include_supersystem=False,
        return_total_data=False,
    )

    # Should compute 1-body in monomer basis, and that's it
    assert r == {
        "all": {1: {((1,), (1,))}},
        "nocp": {1: {((1,), (1,))}},
        "cp": {1: {((1,), (1,))}},
        "vmfc_compute": {1: set()},
        "vmfc_levels": {1: set()},
    }


def test_build_nbody_compute_list_2():
    # Two fragment, only asking for 2 body, nocp
    r = build_nbody_compute_list(
        bsse_type=[BsseEnum.nocp],
        nfragments=2,
        nbodies=[1, 2],
        include_supersystem=False,
        return_total_data=False,
    )

    # Should compute 1-body in monomer basis, dimer in dimer basis
    assert r == {
        "all": {
            1: {((1,), (1,)), ((2,), (2,))},
            2: {((1,), (1,)), ((2,), (2,)), ((1, 2), (1, 2))},
        },
        "nocp": {
            1: {((1,), (1,)), ((2,), (2,))},
            2: {((1,), (1,)), ((2,), (2,)), ((1, 2), (1, 2))},
        },
        "cp": {1: set(), 2: set()},
        "vmfc_compute": {1: set(), 2: set()},
        "vmfc_levels": {1: set(), 2: set()},
    }


def test_build_nbody_compute_list_3():
    # Two fragment, only asking for 2 body, cp correction
    r = build_nbody_compute_list(
        bsse_type=[BsseEnum.cp],
        nfragments=2,
        nbodies=[1, 2],
        include_supersystem=False,
        return_total_data=False,
    )

    # Should compute 1-body in dimer basis, dimer in dimer basis
    assert r == {
        "all": {
            1: {((1,), (1, 2)), ((2,), (1, 2))},
            2: {((1,), (1, 2)), ((2,), (1, 2)), ((1, 2), (1, 2))},
        },
        "nocp": {1: set(), 2: set()},
        "cp": {
            1: {((1,), (1, 2)), ((2,), (1, 2))},
            2: {((1,), (1, 2)), ((2,), (1, 2)), ((1, 2), (1, 2))},
        },
        "vmfc_compute": {1: set(), 2: set()},
        "vmfc_levels": {1: set(), 2: set()},
    }


def test_build_nbody_compute_list_4():
    # Two fragment VMFC
    r = build_nbody_compute_list(
        bsse_type=[BsseEnum.vmfc],
        nfragments=2,
        nbodies=[1, 2],
        include_supersystem=False,
        return_total_data=False,
    )

    assert r == {
        "all": {
            1: {((1,), (1,)), ((2,), (2,))},
            2: {((2,), (1, 2)), ((1,), (1, 2)), ((1, 2), (1, 2))},
        },
        "cp": {1: set(), 2: set()},
        "nocp": {1: set(), 2: set()},
        "vmfc_compute": {
            1: {((1,), (1,)), ((2,), (2,))},
            2: {((2,), (1, 2)), ((1,), (1, 2)), ((1, 2), (1, 2))},
        },
        "vmfc_levels": {
            1: {((1,), (1,)), ((2,), (2,))},
            2: {((2,), (1, 2)), ((1,), (1, 2)), ((1, 2), (1, 2))},
        },
    }
