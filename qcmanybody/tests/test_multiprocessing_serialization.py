"""Test multiprocessing serialization compatibility.

This module identifies and tests serialization issues that prevent multiprocessing
execution in the QCManyBody parallel execution system.
"""

import pickle
from typing import Any, Dict, List, Tuple

import qcelemental as qcel
from qcelemental.models import AtomicInput, Molecule

# Import our local modules
import sys
import os
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)

from qcmanybody import ManyBodyCore, BsseEnum
from qcmanybody.dependency import FragmentDependency, NBodyDependencyGraph
from qcmanybody.parallel import ParallelConfig


class SerializationTester:
    """Test serialization of various QCManyBody objects."""

    @staticmethod
    def create_test_molecule() -> Molecule:
        """Create a simple test molecule with fragments."""
        return qcel.models.Molecule.from_data("""
        O 0.0000000000  0.0000000000 -0.0657755706
        H 0.0000000000 -0.7590619907  0.5219530189
        H 0.0000000000  0.7590619907  0.5219530189
        O 2.5000000000  0.0000000000  0.0000000000
        H 3.2570000000  0.5860000000  0.0000000000
        H 1.7430000000  0.5860000000  0.0000000000
        """, molecular_charge=0, molecular_multiplicity=1, fragments=[[0, 1, 2], [3, 4, 5]])

    @staticmethod
    def test_object_serialization(obj: Any, name: str) -> Tuple[bool, str]:
        """Test if an object can be pickled/unpickled."""
        try:
            # Test pickling
            pickled = pickle.dumps(obj)

            # Test unpickling
            unpickled = pickle.loads(pickled)

            return True, f"✓ {name} serializes successfully"
        except Exception as e:
            return False, f"✗ {name} serialization failed: {e}"

    @staticmethod
    def test_fragment_dependency_serialization():
        """Test FragmentDependency serialization."""
        print("=== Testing FragmentDependency Serialization ===")

        # Create test molecule
        mol = SerializationTester.create_test_molecule()

        # Create FragmentDependency
        try:
            fragment_dep = FragmentDependency("hf", '["hf", [1], [1]]', mol)

            # Test the FragmentDependency object
            success, msg = SerializationTester.test_object_serialization(fragment_dep, "FragmentDependency")
            print(msg)

            if not success:
                # Test individual components
                print("Testing FragmentDependency components:")
                SerializationTester.test_object_serialization(fragment_dep.mc, "  mc (method)")
                SerializationTester.test_object_serialization(fragment_dep.label, "  label")
                SerializationTester.test_object_serialization(fragment_dep.mol, "  mol (Molecule)")

                # Test cached properties
                try:
                    real_atoms = fragment_dep.real_atoms
                    SerializationTester.test_object_serialization(real_atoms, "  real_atoms")
                except:
                    print("  real_atoms: Could not access")

                try:
                    basis_atoms = fragment_dep.basis_atoms
                    SerializationTester.test_object_serialization(basis_atoms, "  basis_atoms")
                except:
                    print("  basis_atoms: Could not access")

        except Exception as e:
            print(f"✗ Could not create FragmentDependency: {e}")

    @staticmethod
    def test_dependency_graph_serialization():
        """Test NBodyDependencyGraph serialization."""
        print("\\n=== Testing NBodyDependencyGraph Serialization ===")

        try:
            # Create a simple compute_map
            mol = SerializationTester.create_test_molecule()
            core = ManyBodyCore(
                molecule=mol,
                bsse_type=[BsseEnum.nocp],
                levels={1: "hf", 2: "hf"},
                return_total_data=False,
                supersystem_ie_only=False,
                embedding_charges={}
            )

            # Get the dependency graph
            dep_graph = core.dependency_graph

            # Test the dependency graph
            success, msg = SerializationTester.test_object_serialization(dep_graph, "NBodyDependencyGraph")
            print(msg)

            if not success:
                # Test components
                print("Testing NBodyDependencyGraph components:")
                try:
                    frag_deps = list(dep_graph.fragment_dependencies)
                    SerializationTester.test_object_serialization(frag_deps, "  fragment_dependencies list")

                    if frag_deps:
                        SerializationTester.test_object_serialization(frag_deps[0], "  first FragmentDependency")

                except Exception as e:
                    print(f"  Could not access fragment_dependencies: {e}")

        except Exception as e:
            print(f"✗ Could not create dependency graph: {e}")

    @staticmethod
    def test_problematic_objects():
        """Test objects that are known to cause serialization issues."""
        print("\\n=== Testing Known Problematic Objects ===")

        # Test dict_keys object
        test_dict = {"a": 1, "b": 2, "c": 3}
        dict_keys_obj = test_dict.keys()

        success, msg = SerializationTester.test_object_serialization(dict_keys_obj, "dict_keys object")
        print(msg)

        # Show the fix - convert to list
        dict_keys_list = list(dict_keys_obj)
        success, msg = SerializationTester.test_object_serialization(dict_keys_list, "dict_keys converted to list")
        print(msg)

    @staticmethod
    def run_comprehensive_serialization_test():
        """Run all serialization tests."""
        print("QCManyBody Multiprocessing Serialization Analysis")
        print("=" * 60)

        SerializationTester.test_fragment_dependency_serialization()
        SerializationTester.test_dependency_graph_serialization()
        SerializationTester.test_problematic_objects()

        print("\\n" + "=" * 60)
        print("Serialization analysis complete.")


if __name__ == "__main__":
    SerializationTester.run_comprehensive_serialization_test()