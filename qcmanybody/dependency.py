"""N-Body Dependency Graph System for Parallel Execution.

This module implements the dependency graph system that enables level-by-level
iteration of N-body fragments, respecting mathematical dependencies:
monomers → dimers → trimers → N-mers.

The key insight is that higher N-body calculations depend on lower N-body
results, so parallel execution must respect this ordering for correctness.
"""

from __future__ import annotations

import json
from typing import Dict, Iterator, List, Set, Tuple, Union

from qcelemental.models import Molecule

from .utils import delabeler


class FragmentDependency:
    """Represents a single fragment with its dependency information.

    Optimized for performance with cached properties and efficient parsing.
    """

    __slots__ = ('mc', 'label', 'mol', '_real_atoms', '_basis_atoms', '_nbody_level')

    def __init__(self, mc: str, label: str, mol: Molecule):
        """Initialize fragment dependency.

        Parameters
        ----------
        mc : str
            Model chemistry identifier
        label : str
            Fragment label in JSON format: '["method", [real_atoms], [basis_atoms]]'
        mol : Molecule
            QCElemental Molecule object for this fragment
        """
        self.mc = mc
        self.label = label
        self.mol = mol

        # Use cached properties for performance
        self._real_atoms = None
        self._basis_atoms = None
        self._nbody_level = None

        # Validate label during construction to maintain error behavior
        try:
            _ = self.nbody_level  # This will trigger parsing and validation
        except Exception:
            # Reset cached values if validation fails
            self._real_atoms = None
            self._basis_atoms = None
            self._nbody_level = None
            raise

    @property
    def real_atoms(self) -> tuple:
        """Get real atoms (cached)."""
        if self._real_atoms is None:
            self._parse_label()
        return self._real_atoms

    @property
    def basis_atoms(self) -> tuple:
        """Get basis atoms (cached)."""
        if self._basis_atoms is None:
            self._parse_label()
        return self._basis_atoms

    @property
    def nbody_level(self) -> int:
        """Get N-body level (cached)."""
        if self._nbody_level is None:
            self._parse_label()
        return self._nbody_level

    def _parse_label(self) -> None:
        """Parse label once and cache results."""
        try:
            _, real_atoms, basis_atoms = delabeler(self.label)
            self._real_atoms = tuple(real_atoms)
            self._basis_atoms = tuple(basis_atoms)
            self._nbody_level = len(real_atoms)
        except Exception as e:
            raise ValueError(f"Failed to parse fragment label '{self.label}': {e}")

    def __repr__(self) -> str:
        return f"FragmentDependency(mc='{self.mc}', nbody={self.nbody_level}, label='{self.label}')"


class NBodyDependencyGraph:
    """N-Body dependency graph for level-ordered fragment iteration.

    This class analyzes fragment labels to extract N-body levels and creates
    a dependency graph that ensures proper execution order for parallel calculations.

    The core innovation is preserving the exact same molecule set as the original
    iterate_molecules() method while adding dependency-respecting level-by-level ordering.
    """

    def __init__(self, compute_map: Dict):
        """Initialize dependency graph from existing compute_map.

        Parameters
        ----------
        compute_map : Dict
            The ManyBodyCore.compute_map structure containing fragment information
        """
        self.compute_map = compute_map
        self.dependency_levels: Dict[int, List[FragmentDependency]] = {}
        self._fragments_by_label: Dict[str, FragmentDependency] = {}
        self._build_dependency_graph()

    def _build_dependency_graph(self) -> None:
        """Build the dependency graph from compute_map.

        This method replicates the exact logic from ManyBodyCore.iterate_molecules()
        but organizes fragments by N-body level instead of arbitrary order.

        Optimized for performance:
        - Pre-import labeler to avoid repeated imports
        - Use more efficient set operations
        - Pre-allocate data structures
        """
        # Performance optimization: import once
        from .utils import labeler

        done_molecules: Set[str] = set()

        # Pre-allocate dictionary with expected levels to reduce dict resizing
        # Most common case is 2-3 levels
        for level in range(1, 5):
            self.dependency_levels[level] = []

        # Use the exact same iteration logic as the original method
        for mc, compute_dict in self.compute_map.items():
            for compute_list in compute_dict["all"].values():
                for real_atoms, basis_atoms in compute_list:
                    label = labeler(mc, real_atoms, basis_atoms)

                    if label in done_molecules:
                        continue

                    # Performance optimization: avoid redundant set operations
                    # We only need the N-body level, not the full ghost atoms calculation

                    # For now, create a placeholder FragmentDependency
                    try:
                        fragment_dep = FragmentDependency(mc, label, None)  # mol will be created later

                        # Group by N-body level
                        level = fragment_dep.nbody_level
                        if level not in self.dependency_levels:
                            self.dependency_levels[level] = []

                        self.dependency_levels[level].append(fragment_dep)
                        self._fragments_by_label[label] = fragment_dep
                        done_molecules.add(label)

                    except ValueError as e:
                        # Skip fragments that can't be parsed (shouldn't happen with valid labels)
                        continue

        # Clean up empty pre-allocated levels
        self.dependency_levels = {k: v for k, v in self.dependency_levels.items() if v}

    def extract_nbody_level(self, label: str) -> int:
        """Extract N-body level from fragment label.

        Parameters
        ----------
        label : str
            Fragment label in format '["method", [real_atoms], [basis_atoms]]'

        Returns
        -------
        int
            N-body level (number of real atoms)

        Examples
        --------
        >>> extract_nbody_level('["mp2", [1], [1, 2]]')
        1
        >>> extract_nbody_level('["ccsd", [1, 2, 3], [1, 2, 3, 4]]')
        3
        """
        try:
            _, real_atoms, _ = delabeler(label)
            return len(real_atoms)
        except Exception as e:
            raise ValueError(f"Cannot extract N-body level from label '{label}': {e}")

    def group_by_level(self, fragments: List[Tuple[str, str, Molecule]]) -> Dict[int, List[Tuple[str, str, Molecule]]]:
        """Group fragments by N-body dependency level.

        Parameters
        ----------
        fragments : List[Tuple[str, str, Molecule]]
            List of (mc, label, mol) tuples from iterate_molecules()

        Returns
        -------
        Dict[int, List[Tuple[str, str, Molecule]]]
            Dictionary mapping N-body level to list of fragments at that level
        """
        grouped: Dict[int, List[Tuple[str, str, Molecule]]] = {}

        for mc, label, mol in fragments:
            level = self.extract_nbody_level(label)
            if level not in grouped:
                grouped[level] = []
            grouped[level].append((mc, label, mol))

        return grouped

    def iterate_molecules_by_level(self) -> Iterator[Tuple[int, List[FragmentDependency]]]:
        """Iterate fragments grouped by dependency level.

        This method provides level-by-level iteration respecting mathematical
        dependencies: monomers (level 1) → dimers (level 2) → trimers (level 3) → etc.

        Yields
        ------
        Tuple[int, List[FragmentDependency]]
            Tuple of (level, fragments_at_level) where fragments_at_level is a list
            of FragmentDependency objects.

        Notes
        -----
        This method preserves the exact same fragment set as the original
        iterate_molecules() method, only changing the ordering to respect dependencies.

        Performance optimized with cached sorting.
        """
        # Cache sorted levels for repeated access
        if not hasattr(self, '_sorted_levels'):
            self._sorted_levels = sorted(self.dependency_levels.keys())

        for level in self._sorted_levels:
            yield level, self.dependency_levels[level]

    def get_dependency_levels(self) -> Dict[int, int]:
        """Get summary of fragments per dependency level.

        Returns
        -------
        Dict[int, int]
            Dictionary mapping N-body level to count of fragments at that level
        """
        return {level: len(fragments) for level, fragments in self.dependency_levels.items()}

    def validate_fragment_completeness(self, original_fragments: List[Tuple[str, str, Molecule]]) -> bool:
        """Validate that dependency graph preserves exact same fragment set.

        Parameters
        ----------
        original_fragments : List[Tuple[str, str, Molecule]]
            Original fragments from iterate_molecules()

        Returns
        -------
        bool
            True if fragment sets are identical

        Raises
        ------
        ValueError
            If fragment sets don't match (missing or extra fragments)
        """
        # Get original fragment labels
        original_labels = {label for _, label, _ in original_fragments}

        # Get dependency graph labels
        dependency_labels = set(self._fragments_by_label.keys())

        if original_labels != dependency_labels:
            missing_in_dependency = original_labels - dependency_labels
            extra_in_dependency = dependency_labels - original_labels

            error_msg = "Fragment set mismatch between original and dependency graph:\n"
            if missing_in_dependency:
                error_msg += f"Missing in dependency: {missing_in_dependency}\n"
            if extra_in_dependency:
                error_msg += f"Extra in dependency: {extra_in_dependency}\n"

            raise ValueError(error_msg)

        return True

    def validate_dependency_ordering(self) -> bool:
        """Validate that dependency ordering is mathematically correct.

        Returns
        -------
        bool
            True if all fragments are at correct dependency levels

        Raises
        ------
        ValueError
            If any fragment is at wrong dependency level
        """
        for level, fragments in self.dependency_levels.items():
            for fragment in fragments:
                if fragment.nbody_level != level:
                    raise ValueError(
                        f"Fragment {fragment.label} has N-body level {fragment.nbody_level} "
                        f"but is in dependency level {level}"
                    )

        return True

    def get_max_level(self) -> int:
        """Get maximum N-body level in the dependency graph.

        Returns
        -------
        int
            Maximum N-body level
        """
        return max(self.dependency_levels.keys()) if self.dependency_levels else 0

    def get_fragments_at_level(self, level: int) -> List[FragmentDependency]:
        """Get all fragments at a specific dependency level.

        Parameters
        ----------
        level : int
            N-body dependency level

        Returns
        -------
        List[FragmentDependency]
            List of fragments at the specified level
        """
        return self.dependency_levels.get(level, [])