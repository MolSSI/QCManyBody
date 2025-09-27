"""Comprehensive tests for N-body dependency graph system.

This test suite validates the dependency graph implementation with ultra-strict
requirements for mathematical correctness and fragment preservation.
"""

import pytest
import numpy as np
from qcelemental.models import Molecule

from qcmanybody.dependency import NBodyDependencyGraph, FragmentDependency
from qcmanybody.utils import labeler, delabeler


class TestFragmentDependency:
    """Test the FragmentDependency class."""

    def test_fragment_dependency_initialization(self):
        """Test FragmentDependency initialization with valid labels."""
        # Test monomer
        label_1b = labeler("hf", (1,), (1,))
        frag = FragmentDependency("hf", label_1b, None)
        assert frag.mc == "hf"
        assert frag.label == label_1b
        assert frag.nbody_level == 1
        assert frag.real_atoms == (1,)
        assert frag.basis_atoms == (1,)

        # Test dimer
        label_2b = labeler("mp2", (1, 2), (1, 2))
        frag = FragmentDependency("mp2", label_2b, None)
        assert frag.nbody_level == 2
        assert frag.real_atoms == (1, 2)
        assert frag.basis_atoms == (1, 2)

        # Test trimer with ghost atom
        label_3b = labeler("ccsd", (1, 2, 3), (1, 2, 3, 4))
        frag = FragmentDependency("ccsd", label_3b, None)
        assert frag.nbody_level == 3
        assert frag.real_atoms == (1, 2, 3)
        assert frag.basis_atoms == (1, 2, 3, 4)

    def test_fragment_dependency_invalid_label(self):
        """Test FragmentDependency with invalid labels."""
        with pytest.raises(ValueError, match="Failed to parse fragment label"):
            FragmentDependency("hf", "invalid_label", None)

    def test_fragment_dependency_repr(self):
        """Test FragmentDependency string representation."""
        label = labeler("mp2", (1, 2), (1, 2, 3))
        frag = FragmentDependency("mp2", label, None)
        repr_str = repr(frag)
        assert "FragmentDependency" in repr_str
        assert "mp2" in repr_str
        assert "nbody=2" in repr_str


class TestNBodyDependencyGraph:
    """Test the NBodyDependencyGraph class."""

    def test_extract_nbody_level(self):
        """Test N-body level extraction from various fragment labels."""
        # Create a minimal dependency graph for testing
        dep_graph = NBodyDependencyGraph({})

        test_cases = [
            (labeler("hf", (1,), (1,)), 1),                    # Monomer
            (labeler("mp2", (1, 2), (1, 2)), 2),               # Dimer
            (labeler("ccsd", (1, 2, 3), (1, 2, 3)), 3),        # Trimer
            (labeler("hf", (2, 4), (1, 2, 3, 4)), 2),          # Dimer with ghosts
            (labeler("mp2", (1, 3, 5), (1, 2, 3, 4, 5)), 3),   # Trimer with ghosts
        ]

        for label, expected_level in test_cases:
            actual_level = dep_graph.extract_nbody_level(label)
            assert actual_level == expected_level, f"Failed for label {label}"

    def test_extract_nbody_level_invalid_label(self):
        """Test N-body level extraction with invalid labels."""
        dep_graph = NBodyDependencyGraph({})

        with pytest.raises(ValueError, match="Cannot extract N-body level"):
            dep_graph.extract_nbody_level("invalid_label")

    def test_group_by_level(self):
        """Test fragment grouping by dependency level."""
        dep_graph = NBodyDependencyGraph({})

        # Create test fragments
        fragments = [
            ("hf", labeler("hf", (1,), (1,)), None),              # Level 1
            ("hf", labeler("hf", (2,), (2,)), None),              # Level 1
            ("mp2", labeler("mp2", (1, 2), (1, 2)), None),        # Level 2
            ("mp2", labeler("mp2", (1, 3), (1, 3)), None),        # Level 2
            ("ccsd", labeler("ccsd", (1, 2, 3), (1, 2, 3)), None), # Level 3
        ]

        grouped = dep_graph.group_by_level(fragments)

        # Verify grouping
        assert len(grouped) == 3  # Three levels: 1, 2, 3
        assert len(grouped[1]) == 2  # Two monomers
        assert len(grouped[2]) == 2  # Two dimers
        assert len(grouped[3]) == 1  # One trimer

        # Verify all fragments preserved
        total_fragments = sum(len(level_frags) for level_frags in grouped.values())
        assert total_fragments == len(fragments)

    def test_dependency_graph_construction_simple(self):
        """Test dependency graph construction with simple compute_map."""
        # Create a simple compute_map mimicking ManyBodyCore structure
        compute_map = {
            "hf": {
                "all": {
                    "1": [((1,), (1,)), ((2,), (2,))],              # Two monomers
                    "2": [((1, 2), (1, 2))],                        # One dimer
                }
            }
        }

        dep_graph = NBodyDependencyGraph(compute_map)

        # Verify dependency levels
        levels = dep_graph.get_dependency_levels()
        assert levels == {1: 2, 2: 1}  # 2 monomers, 1 dimer

        # Verify fragments at each level
        level_1_frags = dep_graph.get_fragments_at_level(1)
        level_2_frags = dep_graph.get_fragments_at_level(2)

        assert len(level_1_frags) == 2
        assert len(level_2_frags) == 1

        # Verify N-body levels are correct
        for frag in level_1_frags:
            assert frag.nbody_level == 1
        for frag in level_2_frags:
            assert frag.nbody_level == 2

    def test_dependency_graph_construction_multilevel(self):
        """Test dependency graph with multi-level calculations."""
        compute_map = {
            "hf": {
                "all": {
                    "1": [((1,), (1,)), ((2,), (2,)), ((3,), (3,))],  # Three monomers
                    "2": [((1, 2), (1, 2)), ((1, 3), (1, 3)), ((2, 3), (2, 3))], # Three dimers
                }
            },
            "mp2": {
                "all": {
                    "3": [((1, 2, 3), (1, 2, 3))],                    # One trimer
                }
            }
        }

        dep_graph = NBodyDependencyGraph(compute_map)

        # Verify dependency levels
        levels = dep_graph.get_dependency_levels()
        assert levels == {1: 3, 2: 3, 3: 1}  # 3 monomers, 3 dimers, 1 trimer

        # Verify max level
        assert dep_graph.get_max_level() == 3

    def test_iterate_molecules_by_level_ordering(self):
        """Test that iterate_molecules_by_level yields in correct dependency order."""
        compute_map = {
            "hf": {
                "all": {
                    "1": [((1,), (1,)), ((2,), (2,))],
                    "2": [((1, 2), (1, 2))],
                    "3": [((1, 2, 3), (1, 2, 3))],
                }
            }
        }

        dep_graph = NBodyDependencyGraph(compute_map)

        # Collect levels in iteration order
        iteration_levels = []
        for level, fragments in dep_graph.iterate_molecules_by_level():
            iteration_levels.append(level)

        # Verify correct dependency order: 1 → 2 → 3
        assert iteration_levels == [1, 2, 3]

    def test_validate_dependency_ordering(self):
        """Test validation of dependency ordering."""
        compute_map = {
            "hf": {
                "all": {
                    "1": [((1,), (1,))],
                    "2": [((1, 2), (1, 2))],
                }
            }
        }

        dep_graph = NBodyDependencyGraph(compute_map)

        # Should pass validation
        assert dep_graph.validate_dependency_ordering() is True

    def test_validate_fragment_completeness_success(self):
        """Test fragment completeness validation with matching sets."""
        compute_map = {
            "hf": {
                "all": {
                    "1": [((1,), (1,)), ((2,), (2,))],
                    "2": [((1, 2), (1, 2))],
                }
            }
        }

        dep_graph = NBodyDependencyGraph(compute_map)

        # Create original fragments that match
        original_fragments = [
            ("hf", labeler("hf", (1,), (1,)), None),
            ("hf", labeler("hf", (2,), (2,)), None),
            ("hf", labeler("hf", (1, 2), (1, 2)), None),
        ]

        # Should pass validation
        assert dep_graph.validate_fragment_completeness(original_fragments) is True

    def test_validate_fragment_completeness_failure(self):
        """Test fragment completeness validation with mismatched sets."""
        compute_map = {
            "hf": {
                "all": {
                    "1": [((1,), (1,))],
                }
            }
        }

        dep_graph = NBodyDependencyGraph(compute_map)

        # Create original fragments that don't match
        original_fragments = [
            ("hf", labeler("hf", (1,), (1,)), None),
            ("hf", labeler("hf", (2,), (2,)), None),  # Extra fragment
        ]

        # Should fail validation
        with pytest.raises(ValueError, match="Fragment set mismatch"):
            dep_graph.validate_fragment_completeness(original_fragments)

    def test_empty_compute_map(self):
        """Test dependency graph with empty compute_map."""
        dep_graph = NBodyDependencyGraph({})

        assert dep_graph.get_dependency_levels() == {}
        assert dep_graph.get_max_level() == 0
        assert list(dep_graph.iterate_molecules_by_level()) == []

    def test_duplicate_fragment_handling(self):
        """Test that duplicate fragments are handled correctly."""
        # This mimics the done_molecules logic in iterate_molecules()
        compute_map = {
            "hf": {
                "all": {
                    "1": [((1,), (1,)), ((1,), (1,))],  # Duplicate fragment
                }
            }
        }

        dep_graph = NBodyDependencyGraph(compute_map)

        # Should only have one fragment (duplicate removed)
        levels = dep_graph.get_dependency_levels()
        assert levels == {1: 1}


class TestDependencyGraphIntegration:
    """Integration tests for dependency graph with real QCManyBody workflows."""

    def test_fragment_preservation_property(self):
        """Property-based test: dependency graph preserves fragment set."""
        # This is a more comprehensive test that would use hypothesis
        # or property-based testing in a full implementation

        # Test with various compute_map structures
        test_cases = [
            # Simple case
            {
                "hf": {
                    "all": {
                        "1": [((1,), (1,))],
                    }
                }
            },
            # Multi-level case
            {
                "hf": {
                    "all": {
                        "1": [((1,), (1,)), ((2,), (2,))],
                        "2": [((1, 2), (1, 2))],
                    }
                }
            },
            # Ghost atoms case
            {
                "mp2": {
                    "all": {
                        "2": [((1, 2), (1, 2, 3))],  # Ghost atom 3
                    }
                }
            }
        ]

        for compute_map in test_cases:
            dep_graph = NBodyDependencyGraph(compute_map)

            # Verify that iteration preserves all fragments
            all_labels_from_iteration = set()
            for level, fragments in dep_graph.iterate_molecules_by_level():
                for fragment_dep in fragments:
                    all_labels_from_iteration.add(fragment_dep.label)

            # Should match internal fragment labels
            assert all_labels_from_iteration == set(dep_graph._fragments_by_label.keys())

    def test_mathematical_dependency_correctness(self):
        """Test that mathematical dependencies are correctly enforced."""
        compute_map = {
            "hf": {
                "all": {
                    "1": [((1,), (1,)), ((2,), (2,)), ((3,), (3,))],
                    "2": [((1, 2), (1, 2)), ((1, 3), (1, 3)), ((2, 3), (2, 3))],
                    "3": [((1, 2, 3), (1, 2, 3))],
                }
            }
        }

        dep_graph = NBodyDependencyGraph(compute_map)

        # Collect all fragments in dependency order
        dependency_order = []
        for level, fragments in dep_graph.iterate_molecules_by_level():
            for fragment_dep in fragments:
                dependency_order.append(fragment_dep.nbody_level)

        # Verify that N-body levels are non-decreasing (dependencies respected)
        for i in range(1, len(dependency_order)):
            assert dependency_order[i] >= dependency_order[i-1], \
                f"Dependency violation: level {dependency_order[i]} before {dependency_order[i-1]}"


@pytest.mark.parametrize("bsse_type", ["cp", "nocp", "vmfc"])
def test_dependency_graph_with_bsse_types(bsse_type):
    """Test dependency graph works with different BSSE treatment types."""
    # This test would be expanded with actual BSSE-specific compute_maps
    # For now, just verify the interface works
    compute_map = {
        "hf": {
            "all": {
                "1": [((1,), (1,))],
                "2": [((1, 2), (1, 2))],
            }
        }
    }

    dep_graph = NBodyDependencyGraph(compute_map)

    # Basic functionality should work regardless of BSSE type
    assert dep_graph.get_max_level() == 2
    assert len(list(dep_graph.iterate_molecules_by_level())) == 2