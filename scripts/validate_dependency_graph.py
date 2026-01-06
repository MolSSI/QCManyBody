#!/usr/bin/env python3
"""Validation script for N-body dependency graph implementation.

This script demonstrates the new dependency graph functionality and validates
that it preserves exact fragment sets while adding dependency ordering.
"""

import sys

# Force fresh import to ensure we get the updated code
for mod in list(sys.modules.keys()):
    if 'qcmanybody' in mod:
        del sys.modules[mod]

from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum


def create_test_molecules():
    """Create test molecules for validation."""

    # Water dimer
    h2o_dimer = Molecule.from_data("""
    units bohr

    O  0.000000  0.000000  0.000000
    H  0.000000  1.431000  1.107000
    H  0.000000 -1.431000  1.107000
    --
    O  0.000000  0.000000  6.000000
    H  0.000000  1.431000  7.107000
    H  0.000000 -1.431000  7.107000
    """)

    # Water trimer
    h2o_trimer = Molecule.from_data("""
    units bohr

    O  0.000000  0.000000  0.000000
    H  0.000000  1.431000  1.107000
    H  0.000000 -1.431000  1.107000
    --
    O  0.000000  0.000000  6.000000
    H  0.000000  1.431000  7.107000
    H  0.000000 -1.431000  7.107000
    --
    O  6.000000  0.000000  3.000000
    H  7.431000  0.000000  4.107000
    H  4.569000  0.000000  4.107000
    """)

    return {"h2o_dimer": h2o_dimer, "h2o_trimer": h2o_trimer}


def validate_fragment_preservation(mbc, system_name):
    """Validate that both iteration methods yield identical fragment sets."""
    print(f"\n=== Validating Fragment Preservation for {system_name} ===")

    # Get fragments from original method
    original_fragments = list(mbc.iterate_molecules())
    original_labels = {label for _, label, _ in original_fragments}

    # Get fragments from new level-ordered method
    level_ordered_fragments = []
    for level, mc, label, mol in mbc.iterate_molecules_by_level():
        level_ordered_fragments.append((mc, label, mol))

    level_ordered_labels = {label for _, label, _ in level_ordered_fragments}

    # Validate
    if original_labels == level_ordered_labels:
        print(f"✅ Fragment sets are identical: {len(original_labels)} fragments")
        return True
    else:
        print(f"❌ Fragment sets differ!")
        print(f"   Original: {len(original_labels)} fragments")
        print(f"   Level-ordered: {len(level_ordered_labels)} fragments")
        print(f"   Missing in level-ordered: {original_labels - level_ordered_labels}")
        print(f"   Extra in level-ordered: {level_ordered_labels - original_labels}")
        return False


def validate_dependency_ordering(mbc, system_name):
    """Validate that dependency ordering is correct."""
    print(f"\n=== Validating Dependency Ordering for {system_name} ===")

    # Collect fragments by level
    fragments_by_level = {}
    dependency_order = []

    for level, mc, label, mol in mbc.iterate_molecules_by_level():
        if level not in fragments_by_level:
            fragments_by_level[level] = []
        fragments_by_level[level].append((mc, label, mol))
        dependency_order.append(level)

    # Check ordering
    is_ordered = True
    for i in range(1, len(dependency_order)):
        if dependency_order[i] < dependency_order[i-1]:
            print(f"❌ Dependency violation: level {dependency_order[i]} after {dependency_order[i-1]}")
            is_ordered = False

    if is_ordered:
        print(f"✅ Dependency ordering is correct")
        print(f"   Levels present: {sorted(fragments_by_level.keys())}")
        for level in sorted(fragments_by_level.keys()):
            print(f"   Level {level}: {len(fragments_by_level[level])} fragments")
        return True
    else:
        print(f"❌ Dependency ordering is incorrect")
        return False


def demonstrate_dependency_graph():
    """Demonstrate dependency graph functionality."""
    print("🚀 QCManyBody N-Body Dependency Graph Validation")
    print("=" * 60)

    molecules = create_test_molecules()
    test_configurations = [
        {
            "name": "Simple 2-body",
            "molecule": molecules["h2o_dimer"],
            "bsse_type": [BsseEnum.cp],
            "levels": {1: "hf", 2: "mp2"},
        },
        {
            "name": "Multi-level 3-body",
            "molecule": molecules["h2o_trimer"],
            "bsse_type": [BsseEnum.cp],
            "levels": {1: "hf", 2: "mp2", 3: "ccsd"},
        },
        {
            "name": "NOCP BSSE",
            "molecule": molecules["h2o_dimer"],
            "bsse_type": [BsseEnum.nocp],
            "levels": {1: "hf", 2: "mp2"},
        }
    ]

    all_passed = True

    for config in test_configurations:
        print(f"\n🧪 Testing Configuration: {config['name']}")
        print("-" * 40)

        try:
            mbc = ManyBodyCore(
                molecule=config["molecule"],
                bsse_type=config["bsse_type"],
                levels=config["levels"],
                return_total_data=True,
                supersystem_ie_only=False,
                embedding_charges={}
            )

            # Validate fragment preservation
            preserve_ok = validate_fragment_preservation(mbc, config["name"])

            # Validate dependency ordering
            order_ok = validate_dependency_ordering(mbc, config["name"])

            # Demonstrate dependency graph properties
            dep_graph = mbc.dependency_graph
            print(f"\n📊 Dependency Graph Properties:")
            print(f"   Max N-body level: {dep_graph.get_max_level()}")
            print(f"   Fragments per level: {dep_graph.get_dependency_levels()}")

            config_passed = preserve_ok and order_ok
            all_passed = all_passed and config_passed

            if config_passed:
                print(f"✅ Configuration '{config['name']}' PASSED")
            else:
                print(f"❌ Configuration '{config['name']}' FAILED")

        except Exception as e:
            print(f"❌ Configuration '{config['name']}' ERROR: {e}")
            all_passed = False

    print(f"\n{'='*60}")
    if all_passed:
        print("🎉 ALL VALIDATION TESTS PASSED!")
        print("✅ N-body dependency graph is working correctly")
        print("✅ Fragment preservation is validated")
        print("✅ Dependency ordering is validated")
        print("✅ Ready for Phase 2 parallel execution implementation")
    else:
        print("❌ SOME VALIDATION TESTS FAILED!")
        print("🔧 Dependency graph implementation needs fixes")

    return all_passed


def demonstrate_new_api():
    """Demonstrate the new API for level-ordered iteration."""
    print(f"\n{'='*60}")
    print("📖 New API Demonstration")
    print("=" * 60)

    molecule = create_test_molecules()["h2o_trimer"]
    mbc = ManyBodyCore(
        molecule=molecule,
        bsse_type=[BsseEnum.cp],
        levels={1: "hf", 2: "mp2", 3: "ccsd"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={}
    )

    print("\n🔄 Original iterate_molecules() - arbitrary order:")
    for i, (mc, label, mol) in enumerate(mbc.iterate_molecules()):
        if i < 5:  # Show first 5
            print(f"   {mc}: {label}")
        elif i == 5:
            print(f"   ... and {len(list(mbc.iterate_molecules())) - 5} more")
            break

    print("\n🎯 New iterate_molecules_by_level() - dependency order:")
    for level, mc, label, mol in mbc.iterate_molecules_by_level():
        print(f"   Level {level}: {mc}: {label}")

    print("\n💡 Key Benefits:")
    print("   ✅ Same molecule set as original method")
    print("   ✅ Level-by-level ordering respects N-body dependencies")
    print("   ✅ Enables safe parallel execution within each level")
    print("   ✅ Backward compatible - original method unchanged")


if __name__ == "__main__":
    success = demonstrate_dependency_graph()
    demonstrate_new_api()

    if success:
        print(f"\n🚀 Phase 1 Task P1-001 Implementation: ✅ COMPLETE")
        print("Ready to proceed with Phase 1 Task P1-002: Integration refinement")
        sys.exit(0)
    else:
        print(f"\n❌ Phase 1 Task P1-001 Implementation: INCOMPLETE")
        print("Fix validation issues before proceeding")
        sys.exit(1)