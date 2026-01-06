#!/usr/bin/env python3
"""
QCManyBody input script generated from HMBE input file: water32_4_by_8_3_4_hmbe.inp

This script sets up a many-body expansion calculation using QCManyBody.
"""

from qcelemental.models import Molecule
from qcmanybody.models import (ManyBodyInput, ManyBodyKeywords,
                               ManyBodySpecification)
from qcmanybody import ManyBodyComputer

def create_molecule():
    """Create the fragmented molecular system."""
    
    # Molecular system with 32 fragments
    symbols = ['O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H']
    
    # Geometry in Bohr (from HMBE input)
    geometry = [2.005979, -4.685521, -1.844031, 1.355527, -4.283326, -1.268374, 1.526601, -4.878857, -2.649665, 2.608892, -2.675449, 2.890175, 2.427245, -3.533869, 2.507619, 2.493361, -2.060559, 2.16575, 0.069765, -3.765074, -0.092294, -0.748399, -4.102246, -0.457193, -0.205709, -3.169585, 0.604654, 0.0, 0.0, 0.0, -0.900139, -0.263684, 0.190917, 0.092177, -0.127148, -0.944228, 2.017875, -5.101032, 1.831973, 2.301866, -5.927774, 1.442003, 1.219745, -4.872587, 1.355497, -2.841906, -3.401469, -1.848014, -2.92869, -2.964138, -1.000996, -2.511626, -2.722254, -2.43607, 3.585351, -2.455495, -1.561401, 2.94275, -3.137891, -1.75538, 4.413623, -2.925805, -1.466523, 2.450667, -1.324316, 0.643786, 2.758818, -1.643552, -0.204364, 1.562815, -1.009682, 0.473627, 1.078891, 0.416333, 5.43754, 1.696153, -0.313147, 5.493008, 1.275802, 0.834104, 4.599135, -2.796958, -0.664606, -0.330683, -3.203523, 0.066005, 0.135301, -2.926777, -0.460567, -1.256827, -3.815253, -2.427004, 2.120311, -4.683296, -2.813738, 2.00556, -3.776899, -1.726033, 1.46962, -1.517278, 1.84158, 2.370512, -2.091013, 1.341235, 2.950781, -1.58883, 1.400826, 1.523844, -1.047178, -2.933996, 2.272746, -0.552524, -2.353539, 2.851206, -1.964328, -2.754123, 2.479411, -4.656323, 1.08899, 0.825096, -4.978444, 1.694752, 0.157623, -5.185808, 1.281009, 1.599048, -3.517569, -1.259768, 4.504909, -4.256747, -0.824819, 4.929953, -3.907261, -1.729748, 3.767692, -0.799077, -1.047125, 4.089732, -0.345323, -0.468214, 4.702266, -1.72862, -0.879239, 4.244624, 4.160764, 2.166572, 2.715175, 4.287608, 3.015355, 2.291257, 3.210581, 2.05839, 2.756142, -1.138153, 4.410088, 0.263741, -0.296342, 4.68125, 0.629877, -1.763663, 4.543124, 0.975969, 3.185867, 3.687994, -0.593526, 3.28441, 4.021231, -1.485418, 2.23966, 3.656005, -0.452463, -0.359419, 3.816258, 3.849684, 0.347573, 4.414256, 3.60722, -0.624062, 3.410303, 3.024218, -3.090156, 3.671481, -1.577749, -2.493968, 3.850117, -0.850509, -3.963404, 3.713109, -1.187962, 1.502374, 1.25824, 2.86289, 1.112074, 0.55707, 2.34111, 1.041895, 2.049705, 2.584017, 0.604867, 4.681388, -1.883436, 0.594428, 5.575712, -2.224471, -0.212207, 4.60275, -1.39106, 3.780276, 0.979797, -0.108064, 3.615404, 0.632822, 0.768666, 3.949462, 1.911673, 0.030538, -2.007336, 1.839099, -5.027128, -2.54765, 2.629216, -5.024783, -1.19034, 2.105807, -4.605689, 3.379245, 2.098633, -4.312302, 3.524446, 2.9781, -3.963471, 3.114449, 1.579866, -3.552701, -2.060018, -1.466237, -3.555814, -1.262512, -0.993527, -3.317574, -1.848445, -1.889886, -4.387672, 1.88149, 0.12892, -6.276787, 2.830236, 0.06493, -6.167178, 1.537069, -0.655066, -5.849038, -3.987328, 0.453972, -2.785206, -3.613701, 1.3271, -2.904706, -3.470783, -0.10907, -3.361742, 0.301391, -0.222093, -2.917628, 1.235922, -0.426827, -2.886593, 0.264685, 0.687556, -3.213298, 3.043168, 0.112457, -2.60588, 3.945524, -0.17369, -2.464093, 2.730091, 0.359177, -1.735628, 0.568276, 2.339609, -3.96175, 1.254414, 2.578278, -4.585032, 0.66538, 2.967531, -3.24585]
    
    # Fragment definitions - each water molecule is a separate fragment
    fragments = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 13, 14], [15, 16, 17], [18, 19, 20], [21, 22, 23], [24, 25, 26], [27, 28, 29], [30, 31, 32], [33, 34, 35], [36, 37, 38], [39, 40, 41], [42, 43, 44], [45, 46, 47], [48, 49, 50], [51, 52, 53], [54, 55, 56], [57, 58, 59], [60, 61, 62], [63, 64, 65], [66, 67, 68], [69, 70, 71], [72, 73, 74], [75, 76, 77], [78, 79, 80], [81, 82, 83], [84, 85, 86], [87, 88, 89], [90, 91, 92], [93, 94, 95]]
    
    # Create the molecule object
    first_frag = {'id': 'PF1', 'fragments': [{'id': 'PF1_SF1', 'symbols': ['O', 'H', 'H'], 'geometry': [2.005979, -4.685521, -1.844031, 1.355527, -4.283326, -1.268374, 1.526601, -4.878857, -2.649665], 'molecular_charge': 0, 'molecular_multiplicity': 1}, {'id': 'PF1_SF2', 'symbols': ['O', 'H', 'H'], 'geometry': [2.608892, -2.675449, 2.890175, 2.427245, -3.533869, 2.507619, 2.493361, -2.060559, 2.16575], 'molecular_charge': 0, 'molecular_multiplicity': 1}, {'id': 'PF1_SF3', 'symbols': ['O', 'H', 'H'], 'geometry': [0.069765, -3.765074, -0.092294, -0.748399, -4.102246, -0.457193, -0.205709, -3.169585, 0.604654], 'molecular_charge': 0, 'molecular_multiplicity': 1}, {'id': 'PF1_SF4', 'symbols': ['O', 'H', 'H'], 'geometry': [0.0, 0.0, 0.0, -0.900139, -0.263684, 0.190917, 0.092177, -0.127148, -0.944228], 'molecular_charge': 0, 'molecular_multiplicity': 1}, {'id': 'PF1_SF5', 'symbols': ['O', 'H', 'H'], 'geometry': [2.017875, -5.101032, 1.831973, 2.301866, -5.927774, 1.442003, 1.219745, -4.872587, 1.355497], 'molecular_charge': 0, 'molecular_multiplicity': 1}, {'id': 'PF1_SF6', 'symbols': ['O', 'H', 'H'], 'geometry': [-2.841906, -3.401469, -1.848014, -2.92869, -2.964138, -1.000996, -2.511626, -2.722254, -2.43607], 'molecular_charge': 0, 'molecular_multiplicity': 1}, {'id': 'PF1_SF7', 'symbols': ['O', 'H', 'H'], 'geometry': [3.585351, -2.455495, -1.561401, 2.94275, -3.137891, -1.75538, 4.413623, -2.925805, -1.466523], 'molecular_charge': 0, 'molecular_multiplicity': 1}, {'id': 'PF1_SF8', 'symbols': ['O', 'H', 'H'], 'geometry': [2.450667, -1.324316, 0.643786, 2.758818, -1.643552, -0.204364, 1.562815, -1.009682, 0.473627], 'molecular_charge': 0, 'molecular_multiplicity': 1}], 'molecular_charge': 0, 'molecular_multiplicity': 1}
    molecule = Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=fragments,
        molecular_charge=first_frag.get("molecular_charge", 0),
        molecular_multiplicity=first_frag.get("molecular_multiplicity", 1)
    )
    
    return molecule


def create_manybody_input():
    """Create the ManyBodyInput object for the calculation."""
    
    molecule = create_molecule()
    
    # Many-body keywords
    keywords = ManyBodyKeywords(
        bsse_type=["cp"],  # Counterpoise correction
        max_nbody=4,  # Maximum n-body level from HMBE
        return_total_data=True,  # Return total energies/properties
        supersystem_ie_only=False  # Compute all n-body levels
    )
    
    # QC specification for the quantum chemistry method
    specification = ManyBodySpecification(
        driver="energy",  # energy, gradient, or hessian
        keywords=keywords,
        specification={
            "model": {
                "method": "mp2",
                "basis": "aug-cc-pVDZ"
            },
            "driver": "energy",
            "program": "psi4",
            "keywords": {'scf_type': 'df', 'maxiter': 100, 'e_convergence': 1e-06, 'd_convergence': 1e-06},
            "protocols": {
                "stdout": True
            }
        }
    )
    
    # Create the full input
    manybody_input = ManyBodyInput(
        molecule=molecule,
        specification=specification
    )
    
    return manybody_input


def run_calculation():
    """Run the many-body expansion calculation."""
    
    # Create the input
    manybody_input = create_manybody_input()
    
    # Create the computer and run the calculation
    computer = ManyBodyComputer.from_manybodyinput(manybody_input)
    
    # The computer object can now be used to:
    # 1. Generate individual QC calculations: computer.plan()
    # 2. Run calculations through QCEngine: computer.compute()
    # 3. Get results: computer.get_results(...)
    
    print("ManyBody calculation setup complete!")
    n_atoms = len(computer.molecule.symbols)
    n_frags = computer.nfragments
    print(f"System: {n_atoms} atoms in {n_frags} fragments")
    print(f"Max n-body level: {computer.max_nbody}")
    print(f"BSSE treatment: {computer.bsse_type}")
    spec = manybody_input.specification.specification
    print(f"Method: {spec['model']['method']}")
    print(f"Basis: {spec['model']['basis']}")
    
    return computer


if __name__ == "__main__":
    # Run the setup
    computer = run_calculation()
    
    # Optionally, you can inspect the calculation plan:
    print("\nCalculation plan:")
    plan_text, plan_dict = computer.qcmb_core.format_calc_plan()
    print(plan_text)
    
    # Example of running QC calculations (requires QCEngine):
    # computer.compute()  # This would run all QC calculations
    # results = computer.get_results(...)  # Get the results
