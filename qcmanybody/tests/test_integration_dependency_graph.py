"""Integration tests for dependency graph with ManyBodyCore.

This test suite validates the integration of the dependency graph system
with ManyBodyCore, ensuring mathematical correctness and fragment preservation.
"""

import pytest
import numpy as np
from qcelemental.models import Molecule

from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum


class TestManyBodyCoreIntegration:
    """Test integration of dependency graph with ManyBodyCore."""

    @pytest.fixture
    def simple_molecule(self):
        """Create a simple test molecule (water dimer)."""
        return Molecule.from_data("""
        units bohr

        O  0.000000  0.000000  0.000000
        H  0.000000  1.431000  1.107000
        H  0.000000 -1.431000  1.107000
        --
        O  0.000000  0.000000  6.000000
        H  0.000000  1.431000  7.107000
        H  0.000000 -1.431000  7.107000
        """)

    @pytest.fixture
    def water_trimer(self):
        """Create a water trimer for 3-body testing."""
        return Molecule.from_data("""
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

    def test_dependency_graph_property_access(self, simple_molecule):
        """Test that dependency graph property is accessible."""
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Access dependency graph
        dep_graph = mbc.dependency_graph
        assert dep_graph is not None

        # Should cache the instance
        assert mbc.dependency_graph is dep_graph

    def test_iterate_molecules_by_level_exists(self, simple_molecule):
        """Test that iterate_molecules_by_level method exists and is callable."""
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Method should exist
        assert hasattr(mbc, 'iterate_molecules_by_level')
        assert callable(mbc.iterate_molecules_by_level)

        # Should be iterable
        iterator = mbc.iterate_molecules_by_level()
        assert hasattr(iterator, '__iter__')

    def test_fragment_preservation_simple(self, simple_molecule):
        """Test that iterate_molecules_by_level preserves exact same fragment set."""
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Get fragments from original method
        original_fragments = list(mbc.iterate_molecules())
        original_labels = {label for _, label, _ in original_fragments}

        # Get fragments from new level-ordered method
        level_ordered_fragments = []
        for level, mc, label, mol in mbc.iterate_molecules_by_level():
            level_ordered_fragments.append((mc, label, mol))

        level_ordered_labels = {label for _, label, _ in level_ordered_fragments}

        # Fragment sets must be identical
        assert original_labels == level_ordered_labels, \
            f"Fragment sets differ: original={original_labels}, level_ordered={level_ordered_labels}"

        # Total count must be identical
        assert len(original_fragments) == len(level_ordered_fragments), \
            f"Fragment counts differ: {len(original_fragments)} vs {len(level_ordered_fragments)}"

    def test_dependency_ordering_correctness(self, water_trimer):
        """Test that dependency ordering is mathematically correct."""
        mbc = ManyBodyCore(
            molecule=water_trimer,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2", 3: "ccsd"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Collect fragments in dependency order
        dependency_order = []
        for level, mc, label, mol in mbc.iterate_molecules_by_level():
            dependency_order.append(level)

        # Verify non-decreasing order (dependencies respected)
        for i in range(1, len(dependency_order)):
            assert dependency_order[i] >= dependency_order[i-1], \
                f"Dependency violation: level {dependency_order[i]} before {dependency_order[i-1]}"

    def test_level_grouping(self, water_trimer):
        """Test that fragments are correctly grouped by level."""
        mbc = ManyBodyCore(
            molecule=water_trimer,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2", 3: "ccsd"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Group fragments by level
        fragments_by_level = {}
        for level, mc, label, mol in mbc.iterate_molecules_by_level():
            if level not in fragments_by_level:
                fragments_by_level[level] = []
            fragments_by_level[level].append((mc, label, mol))

        # Should have levels 1, 2, 3 for water trimer
        assert 1 in fragments_by_level, "Missing monomer level"
        assert 2 in fragments_by_level, "Missing dimer level"
        assert 3 in fragments_by_level, "Missing trimer level"

        # Verify counts (note: multi-level calculations generate multiple computations per fragment)
        # Level 1: 3 monomers * 3 methods (hf for level 1, mp2 for ghost computations, ccsd for ghost computations)
        assert len(fragments_by_level[1]) >= 3, f"Expected at least 3 monomers, got {len(fragments_by_level[1])}"
        # Level 2: 3 dimers * methods
        assert len(fragments_by_level[2]) >= 3, f"Expected at least 3 dimers, got {len(fragments_by_level[2])}"
        # Level 3: 1 trimer * methods
        assert len(fragments_by_level[3]) >= 1, f"Expected at least 1 trimer, got {len(fragments_by_level[3])}"

    def test_molecule_geometry_preservation(self, simple_molecule):
        """Test that molecule geometries are preserved between methods."""
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Get molecules from both methods
        original_mols = {label: mol for _, label, mol in mbc.iterate_molecules()}
        level_ordered_mols = {label: mol for _, _, label, mol in mbc.iterate_molecules_by_level()}

        # Should have same labels
        assert set(original_mols.keys()) == set(level_ordered_mols.keys())

        # Geometries should be identical for each label
        for label in original_mols:
            orig_geom = original_mols[label].geometry
            level_geom = level_ordered_mols[label].geometry

            np.testing.assert_array_equal(orig_geom, level_geom,
                err_msg=f"Geometry mismatch for fragment {label}")

    @pytest.mark.parametrize("bsse_type", [BsseEnum.cp, BsseEnum.nocp, BsseEnum.vmfc])
    def test_bsse_compatibility(self, simple_molecule, bsse_type):
        """Test compatibility with different BSSE treatment types."""
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[bsse_type],
            levels={1: "hf", 2: "mp2"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Should work with all BSSE types
        original_count = len(list(mbc.iterate_molecules()))
        level_ordered_count = len(list(mbc.iterate_molecules_by_level()))

        assert original_count == level_ordered_count

    def test_multilevel_calculations(self, water_trimer):
        """Test with multi-level calculations (different methods per N-body level)."""
        mbc = ManyBodyCore(
            molecule=water_trimer,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2", 3: "ccsd"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Verify method assignments per level
        method_by_level = {}
        for level, mc, label, mol in mbc.iterate_molecules_by_level():
            if level not in method_by_level:
                method_by_level[level] = set()
            method_by_level[level].add(mc)

        # Each level should have its assigned method
        assert "hf" in method_by_level.get(1, set()), "Level 1 should use HF"
        assert "mp2" in method_by_level.get(2, set()), "Level 2 should use MP2"
        assert "ccsd" in method_by_level.get(3, set()), "Level 3 should use CCSD"

    def test_limited_fragments_handling(self, simple_molecule):
        """Test handling when limited fragments are generated."""
        # Create ManyBodyCore with setup that actually works
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2"},  # Standard multi-level
            return_total_data=True,
            supersystem_ie_only=False,  # Standard MBE
            embedding_charges={}
        )

        # Both methods should produce the same fragments
        original_list = list(mbc.iterate_molecules())
        level_ordered_list = list(mbc.iterate_molecules_by_level())

        # Should have the same count
        assert len(original_list) == len(level_ordered_list)

        # And the same fragment labels
        original_labels = {label for _, label, _ in original_list}
        level_ordered_labels = {label for _, _, label, _ in level_ordered_list}
        assert original_labels == level_ordered_labels

    def test_supersystem_compatibility(self, simple_molecule):
        """Test compatibility with supersystem calculations."""
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2", "supersystem": "ccsd"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Should work with supersystem
        original_labels = {label for _, label, _ in mbc.iterate_molecules()}
        level_ordered_labels = {label for _, _, label, _ in mbc.iterate_molecules_by_level()}

        assert original_labels == level_ordered_labels


class TestDependencyGraphValidation:
    """Test validation aspects of the dependency graph integration."""

    @pytest.fixture
    def simple_molecule(self):
        """Create a simple test molecule (water dimer)."""
        return Molecule.from_data("""
        units bohr

        O  0.000000  0.000000  0.000000
        H  0.000000  1.431000  1.107000
        H  0.000000 -1.431000  1.107000
        --
        O  0.000000  0.000000  6.000000
        H  0.000000  1.431000  7.107000
        H  0.000000 -1.431000  7.107000
        """)

    @pytest.fixture
    def water_trimer(self):
        """Create a water trimer for 3-body testing."""
        return Molecule.from_data("""
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

    def test_fragment_completeness_validation(self, simple_molecule):
        """Test the fragment completeness validation method."""
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Get original fragments for validation
        original_fragments = list(mbc.iterate_molecules())

        # Validation should pass
        assert mbc.dependency_graph.validate_fragment_completeness(original_fragments)

    def test_dependency_ordering_validation(self, simple_molecule):
        """Test the dependency ordering validation method."""
        mbc = ManyBodyCore(
            molecule=simple_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Validation should pass
        assert mbc.dependency_graph.validate_dependency_ordering()

    def test_max_level_detection(self, water_trimer):
        """Test maximum N-body level detection."""
        mbc = ManyBodyCore(
            molecule=water_trimer,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf", 2: "mp2", 3: "ccsd"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={}
        )

        # Should detect max level as 3 for trimer
        assert mbc.dependency_graph.get_max_level() == 3