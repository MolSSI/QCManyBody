from __future__ import annotations

import logging
import math
from typing import Set, Dict, Tuple, Iterable, Union, Literal, Mapping, Any

import numpy as np
from qcelemental.models import DriverEnum, Molecule

from qcmanybody.assemble import _sum_cluster_ptype_data
from qcmanybody.builder import build_nbody_compute_list
from qcmanybody.models import BsseEnum, FragBasIndex
from qcmanybody.utils import delabeler, labeler, print_nbody_energy

logger = logging.getLogger(__name__)


class ManybodyCalculator:
    def __init__(
        self,
        molecule: Molecule,
        bsse_type: Iterable[BsseEnum],
        levels: Mapping[Union[int, Literal["supersystem"]], str],
        return_total_data: bool,
    ):

        # TODO
        self.embedding_charges = {}
        self.driver = DriverEnum.energy

        self.molecule = molecule
        self.bsse_type = set(bsse_type)
        self.return_total_data = return_total_data
        self.nfragments = len(molecule.fragments)

        self.levels = levels

        # Levels without supersystem
        self.levels_no_ss = {int(k): v for k, v in levels.items() if k != "supersystem"}

        # Just a set of all the modelchems
        self.mc_levels = set(self.levels.values())

        self.max_nbody = max(self.levels_no_ss.keys())

        if len(self.bsse_type) == 0:
            raise ValueError("No BSSE correction specified")

        if BsseEnum.vmfc in self.bsse_type:
            self.return_bsse_type = "vmfc"
        elif BsseEnum.cp in self.bsse_type:
            self.return_bsse_type = "cp"
        elif BsseEnum.nocp in self.bsse_type:
            self.return_bsse_type = "nocp"
        else:
            raise RuntimeError("Cannot figure out return_bsse_type. Please post this error on github.")

        ###############################
        # Build nbodies_per_mc_level
        # TODO - use Lori's code
        max_level = max(self.levels_no_ss.keys())

        if set(range(1, max_level + 1)) != set(self.levels_no_ss.keys()):
            raise ValueError(f"Levels must be contiguous from 1 to {max_level}")

        self.nbodies_per_mc_level: Dict[str, list] = {
            mc_level: [] for mc_level in self.mc_levels
        }
        for k, v in self.levels_no_ss.items():
            self.nbodies_per_mc_level[v].append(k)

        self.nbodies_per_mc_level = {
            k: sorted(v) for k, v in self.nbodies_per_mc_level.items()
        }

        # Supersystem is always at the end
        if "supersystem" in levels:
            ss_mc = levels["supersystem"]
            self.nbodies_per_mc_level[ss_mc].append("supersystem")

        # To be built on the fly
        self.mc_compute_dict = None

    @property
    def has_supersystem(self) -> bool:
        return "supersystem" in self.levels

    @property
    def compute_map(self) -> Dict[str, Dict[str, Dict[int, Set[FragBasIndex]]]]:

        if self.mc_compute_dict is not None:
            return self.mc_compute_dict

        # Build the compute lists
        self.mc_compute_dict = {}

        for mc in self.mc_levels:
            nbodies = self.nbodies_per_mc_level[mc]
            self.mc_compute_dict[mc] = build_nbody_compute_list(
                self.bsse_type,
                self.nfragments,
                nbodies,
                self.return_total_data,
                self.max_nbody,
            )

        return self.mc_compute_dict

    def iterate_molecules(self, orient: bool = True) -> Tuple[str, str, Molecule]:
        """Iterate over all the molecules needed for the computation.

        Yields model chemistry, label, and molecule.
        """

        done_molecules = set()

        for mc, compute_dict in self.compute_map.items():

            # TODO - this is a bit of a hack. Lots of duplication when reaching higher nbody
            for compute_list in compute_dict["all"].values():
                for real_atoms, basis_atoms in compute_list:
                    label = labeler(mc, real_atoms, basis_atoms)
                    if label in done_molecules:
                        continue

                    ghost_atoms = list(set(basis_atoms) - set(real_atoms))

                    # Shift to zero-indexing
                    real_atoms_0 = [x - 1 for x in real_atoms]
                    ghost_atoms_0 = [x - 1 for x in ghost_atoms]
                    mol = self.molecule.get_fragment(
                        real_atoms_0, ghost_atoms_0, orient=orient, group_fragments=True
                    )

                    # if self.embedding_charges:
                    #    embedding_frags = list(set(range(1, self.nfragments + 1)) - set(pair[1]))
                    #    charges = []
                    #    for frag in embedding_frags:
                    #        positions = self.molecule.extract_subsets(frag).geometry().np.tolist()
                    #        charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[frag])])
                    #    data['keywords']['function_kwargs'].update({'external_potentials': charges})

                    done_molecules.add(label)
                    yield mc, label, mol

    def _assemble_nbody_components(
        self,
        ptype: DriverEnum,
        component_results: Dict[str, Union[float, np.ndarray]],
    ) -> Dict[str, Any]:
        """Assembles N-body components for a single derivative level and a single model chemistry level
        into interaction quantities according to requested BSSE treatment(s).
        """

        # which level are we assembling?
        delabled = [delabeler(k) for k in component_results.keys()]
        mc_level_labels = {x[0] for x in delabled}

        if len(mc_level_labels) != 1:
            raise RuntimeError(
                f"Multiple model chemistries passed into _assemble_nbody_components: {mc_level_labels}"
            )

        mc_level = mc_level_labels.pop()
        if mc_level not in self.mc_levels:
            raise RuntimeError(
                f"Model chemistry {mc_level} not found in {self.mc_levels}"
            )

        # get the range of nbodies and the required calculations for this level
        nbodies = self.nbodies_per_mc_level[mc_level]
        max_nbody = nbodies[-1] # TODO - or max()?
        compute_dict = self.compute_map[mc_level]

        # Build size and slices dictionaries
        fragment_size_dict = {}
        fragment_slice_dict = {}
        iat = 0
        for ifr in range(1, self.nfragments + 1):
            nat = len(self.molecule.fragments[ifr - 1])
            fragment_size_dict[ifr] = nat
            fragment_slice_dict[ifr] = slice(iat, iat + nat)
            iat += nat

        # Final dictionaries
        if ptype == DriverEnum.energy:
            cp_by_level = {n: 0.0 for n in range(1, nbodies[-1] + 1)}
            nocp_by_level = {n: 0.0 for n in range(1, nbodies[-1] + 1)}
            vmfc_by_level = {n: 0.0 for n in range(1, nbodies[-1] + 1)}

            cp_body_dict = {n: 0.0 for n in range(1, nbodies[-1] + 1)}
            nocp_body_dict = {n: 0.0 for n in range(1, nbodies[-1] + 1)}
            vmfc_body_dict = {n: 0.0 for n in range(1, nbodies[-1] + 1)}

        else:
            nat = sum(fragment_size_dict.values())
            if ptype == DriverEnum.gradient:
                arr_shape = (nat, 3)
            elif ptype == DriverEnum.hessian:
                arr_shape = (nat * 3, nat * 3)
            else:
                raise ValueError(f"Invalid ptype {ptype}")

            cp_by_level = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}
            nocp_by_level = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}
            vmfc_by_level = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}

            cp_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}
            nocp_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}
            vmfc_body_dict = {n: np.zeros(arr_shape) for n in range(1, nbodies[-1] + 1)}

        # Sum up all of the levels
        # * compute_dict[bt][nb] holds all the computations needed to compute nb
        #   *not* all the nb-level computations, so build the latter
        cp_compute_list = {nb: set() for nb in range(1, nbodies[-1] + 1)}
        nocp_compute_list = {nb: set() for nb in range(1, nbodies[-1] + 1)}

        for nb in nbodies:
            for v in compute_dict["cp"][nb]:
                if len(v[1]) != 1:
                    cp_compute_list[len(v[0])].add(v)
            for w in compute_dict["nocp"][nb]:
                nocp_compute_list[len(w[0])].add(w)

        for nb in range(1, nbodies[-1] + 1):
            cp_by_level[nb] = _sum_cluster_ptype_data(
                ptype,
                component_results,
                cp_compute_list[nb],
                fragment_slice_dict,
                fragment_size_dict,
                mc_level,
            )
            nocp_by_level[nb] = _sum_cluster_ptype_data(
                ptype,
                component_results,
                nocp_compute_list[nb],
                fragment_slice_dict,
                fragment_size_dict,
                mc_level,
            )
            if nb in compute_dict["vmfc_levels"]:
                vmfc_by_level[nb] = _sum_cluster_ptype_data(
                    ptype,
                    component_results,
                    compute_dict["vmfc_levels"][nb],
                    fragment_slice_dict,
                    fragment_size_dict,
                    mc_level,
                    vmfc=True,
                    nb=nb,
                )

        # Extract data for monomers in monomer basis for CP total data
        if 1 in nbodies:
            monomers_in_monomer_basis = [
                v for v in compute_dict["nocp"][1] if len(v[1]) == 1
            ]

            if ptype == "energy":
                monomer_energy_list = [
                    component_results[labeler(mc_level, m[0], m[1])]
                    for m in monomers_in_monomer_basis
                ]
                monomer_sum = sum(monomer_energy_list)
            else:
                monomer_sum = _sum_cluster_ptype_data(
                    ptype,
                    component_results,
                    set(monomers_in_monomer_basis),
                    fragment_slice_dict,
                    fragment_size_dict,
                    mc_level,
                )
        else:
            # TODO - fixme
            monomer_sum = 0.0

        nbody_dict = {}

        # Compute cp
        if "cp" in self.bsse_type:
            for nb in range(1, nbodies[-1] + 1):
                if nb == self.nfragments:
                    if ptype == "energy":
                        cp_body_dict[nb] = cp_by_level[nb] - bsse
                    else:
                        cp_body_dict[nb][:] = cp_by_level[nb] - bsse
                    continue

                for k in range(1, nb + 1):
                    take_nk = math.comb(self.nfragments - k - 1, nb - k)
                    sign = (-1) ** (nb - k)
                    cp_body_dict[nb] += take_nk * sign * cp_by_level[k]

                if nb == 1:
                    bsse = cp_body_dict[nb] - monomer_sum
                    if ptype == "energy":
                        cp_body_dict[nb] = monomer_sum
                    else:
                        cp_body_dict[nb] = monomer_sum.copy()
                else:
                    cp_body_dict[nb] -= bsse

            if ptype == "energy":
                print_nbody_energy(
                    cp_body_dict,
                    "Counterpoise Corrected (CP)",
                    self.nfragments,
                    self.embedding_charges,
                )

                if monomer_sum != 0.0:
                    nbody_dict["CP-CORRECTED TOTAL ENERGY"] = cp_body_dict[max_nbody]
                nbody_dict["CP-CORRECTED INTERACTION ENERGY"] = cp_body_dict[max_nbody] - cp_body_dict[1]

                for nb in nbodies[1:]:
                    nbody_dict[f"CP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY"] = cp_body_dict[nb] - cp_body_dict[1]
                    nbody_dict[f"CP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY"] = cp_body_dict[nb] - cp_body_dict[nb - 1]
                for nb in nbodies:
                    nbody_dict[f"CP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY"] = cp_body_dict[nb]

        # Compute nocp
        if "nocp" in self.bsse_type:
            for nb in range(1, nbodies[-1] + 1):
                if nb == self.nfragments:
                    if ptype == "energy":
                        nocp_body_dict[nb] = nocp_by_level[nb]
                    else:
                        nocp_body_dict[nb][:] = nocp_by_level[nb]
                    continue

                for k in range(1, nb + 1):
                    take_nk = math.comb(self.nfragments - k - 1, nb - k)
                    sign = (-1) ** (nb - k)
                    nocp_body_dict[nb] += take_nk * sign * nocp_by_level[k]

            if ptype == "energy":
                print_nbody_energy(
                    nocp_body_dict,
                    "Non-Counterpoise Corrected (NoCP)",
                    self.nfragments,
                    self.embedding_charges,
                )

                nbody_dict["NOCP-CORRECTED TOTAL ENERGY"] = nocp_body_dict[max_nbody]
                nbody_dict["NOCP-CORRECTED INTERACTION ENERGY"] = nocp_body_dict[max_nbody] - nocp_body_dict[1]

                for nb in nbodies[1:]:
                    nbody_dict[
                        f"NOCP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY"
                    ] = (nocp_body_dict[nb] - nocp_body_dict[1])
                    nbody_dict[f"NOCP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY"] = (
                        nocp_body_dict[nb] - nocp_body_dict[nb - 1]
                    )
                for nb in nbodies:
                    nbody_dict[f"NOCP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY"] = (
                        nocp_body_dict[nb]
                    )

        # Compute vmfc
        if "vmfc" in self.bsse_type:
            for nb in nbodies:
                if ptype == "energy":
                    for k in range(1, nb + 1):
                        vmfc_body_dict[nb] += vmfc_by_level[k]

                else:
                    if nb > 1:
                        vmfc_body_dict[nb] = vmfc_by_level[nb - 1]
                    vmfc_body_dict[nb] += vmfc_by_level[nb]

            if ptype == "energy":
                print_nbody_energy(
                    vmfc_body_dict,
                    "Valiron-Mayer Function Counterpoise (VMFC)",
                    self.nfragments,
                    self.embedding_charges,
                )

                vmfc_interaction_energy = (
                    vmfc_body_dict[max_nbody] - vmfc_body_dict[1]
                )
                nbody_dict["VMFC-CORRECTED TOTAL ENERGY"] = vmfc_body_dict[
                    max_nbody
                ]
                nbody_dict["VMFC-CORRECTED INTERACTION ENERGY"] = (
                    vmfc_interaction_energy
                )

                for nb in nbodies[1:]:
                    nbody_dict[
                        f"VMFC-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY"
                    ] = (vmfc_body_dict[nb] - vmfc_body_dict[1])
                    nbody_dict[f"VMFC-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY"] = (
                        vmfc_body_dict[nb] - vmfc_body_dict[nb - 1]
                    )
                for nb in nbodies:
                    nbody_dict[f"VMFC-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY"] = (
                        vmfc_body_dict[nb]
                    )

        # Collect specific and generalized returns
        results = {
            f"cp_{ptype}_body_dict": {f"{nb}cp": j for nb, j in cp_body_dict.items()},
            f"nocp_{ptype}_body_dict": {
                f"{nb}nocp": j for nb, j in nocp_body_dict.items()
            },
            f"vmfc_{ptype}_body_dict": {
                f"{nb}vmfc": j for nb, j in vmfc_body_dict.items()
            },
        }

        if ptype == "energy":
            results["nbody"] = nbody_dict

        # TODO - not well defined?
        #return_bsse_type = self.bsse_type[0]


        if self.return_bsse_type == "cp":
            results[f"{ptype}_body_dict"] = cp_body_dict
        elif self.return_bsse_type == "nocp":
            results[f"{ptype}_body_dict"] = nocp_body_dict
        elif self.return_bsse_type == "vmfc":
            results[f"{ptype}_body_dict"] = vmfc_body_dict
        else:
            raise RuntimeError(
                "N-Body Wrapper: Invalid return type. Should never be here, please post this error on github."
            )

        if ptype == "energy":
            piece = results[f"{ptype}_body_dict"][max_nbody]
        else:
            piece = results[f"{ptype}_body_dict"][max_nbody].copy()

        if self.return_total_data:
            results[f"ret_{ptype}"] = piece
        else:
            results[f"ret_{ptype}"] = piece
            results[f"ret_{ptype}"] -= results[f"{ptype}_body_dict"][1]

        results["ret_ptype"] = results[f"ret_{ptype}"]

        return results

    def _analyze(
            self,
            ptype: DriverEnum,
            component_results: Dict[str, Union[float, np.ndarray]],
    ):

        natoms = len(self.molecule.symbols)

        # Initialize with zeros
        energy_result, gradient_result, hessian_result = 0, None, None
        energy_body_contribution = {b: {} for b in self.bsse_type}
        energy_body_dict = {b: {} for b in self.bsse_type}
        if ptype in [DriverEnum.gradient, DriverEnum.hessian]:
            gradient_result = np.zeros((natoms, 3))
        if ptype == DriverEnum.hessian:
            hessian_result = np.zeros((natoms * 3, natoms * 3))

        # Order the levels by highest nbody (??)
        # TODO - right?
        levels = sorted(self.levels_no_ss.items(), reverse=True)

        for nb, label in levels:

            nbody_list2 = self.nbodies_per_mc_level[label]

            cresults = {k: v for k, v in component_results.items() if delabeler(k)[0] == label}
            results = self._assemble_nbody_components(ptype, cresults)

            for n in nbody_list2[::-1]:
                energy_bsse_dict = {b: 0 for b in self.bsse_type}

                for m in range(n - 1, n + 1):
                    if m == 0:
                        continue
                    # Subtract the (n-1)-body contribution from the n-body contribution to get the n-body effect
                    sign = (-1) ** (1 - m // n)
                    for b in self.bsse_type:
                        energy_bsse_dict[b] += (
                                sign
                                * results["%s_energy_body_dict" % b.lower()][
                                    "%i%s" % (m, b.lower())
                                    ]
                        )

                    if ptype == "hessian":
                        hessian_result += sign * results[f"{ptype}_body_dict"][m]
                        gradient_result += sign * results["gradient_body_dict"][m]
                        if n == 1:
                            hessian1 = results[f"{ptype}_body_dict"][n]
                            gradient1 = results["gradient_body_dict"][n]

                    elif ptype == "gradient":
                        gradient_result += sign * results[f"{ptype}_body_dict"][m]
                        # Keep 1-body contribution to compute interaction data
                        if n == 1:
                            gradient1 = results[f"{ptype}_body_dict"][n]

                energy_result += energy_bsse_dict[self.return_bsse_type]
                for b in self.bsse_type:
                    energy_body_contribution[b][n] = energy_bsse_dict[b]

        if self.has_supersystem:
            # Get the MC label for supersystem tasks
            supersystem_mc_level = self.levels.get('supersystem', None)

            # Super system recovers higher order effects at a lower level
            frag_range = tuple(range(1, self.nfragments + 1))

            supersystem_result = component_results[labeler(supersystem_mc_level, frag_range, frag_range)]

            max_nbody = max(self.levels_no_ss.keys())
            print("FOUND SUPERSYSTEM", supersystem_result)

            # Compute components at supersystem level of theory
            self.nbodies_per_mc_level.append(levels)
            component_result = {
                k: v for k, v in self.task_list.items() if k.startswith(str(sup_level))
            }
            components = self.prepare_results(results=component_result, client=client)

            energy_result += (
                    supersystem_result.properties.return_energy
                    - components["energy_body_dict"][self.max_nbody]
            )
            for b in self.bsse_type:
                energy_body_contribution[b][self.molecule.nfragments()] = (
                        supersystem_result.properties.return_energy
                        - components["energy_body_dict"][self.max_nbody]
                )

            if ptype == "hessian":
                gradient_result += (
                        supersystem_result.extras.qcvars["CURRENT GRADIENT"]
                        - components["gradient_body_dict"][self.max_nbody]
                )
                hessian_result += (
                        supersystem_result.return_result
                        - components[f"{ptype}_body_dict"][self.max_nbody]
                )

            elif ptype == "gradient":
                gradient_result += (
                        np.array(supersystem_result.return_result).reshape((-1, 3))
                        - components[f"{ptype}_body_dict"][self.max_nbody]
                )

        for b in self.bsse_type:
            for n in energy_body_contribution[b]:
                energy_body_dict[b][n] = sum(
                    [
                        energy_body_contribution[b][i]
                        for i in range(1, n + 1)
                        if i in energy_body_contribution[b]
                    ]
                )

        is_embedded = self.embedding_charges
        for b in self.bsse_type:
            print_nbody_energy(
                energy_body_dict[b],
                f"{b.upper()}-corrected multilevel many-body expansion",
                self.nfragments,
                is_embedded,
            )

        if not self.return_total_data:
            # Remove monomer cotribution for interaction data
            energy_result -= energy_body_dict[self.return_bsse_type][1]
            if ptype in ["gradient", "hessian"]:
                gradient_result -= gradient1
            if ptype == "hessian":
                hessian_result -= hessian1

        energy_body_dict = {
            str(k) + b: v for b in energy_body_dict for k, v in energy_body_dict[b].items()
        }

        nbody_results = {
            "ret_energy": energy_result,
            "ret_ptype": locals()[ptype + "_result"],
            "energy_body_dict": energy_body_dict,
        }
        return nbody_results

    def analyze(
            self,
            component_results: Dict[str, Dict[str, Union[float, np.ndarray]]], # component_results[label]['energy'] = 1.23
    ):

        # reorganize to component_results_inv['energy'][label] = 1.23
        component_results_inv = {'energy': {}, 'gradient': {}, 'hessian': {}}
        for label, data in component_results.items():
                component_results_inv['energy'][label] = data.get('energy', None)
                component_results_inv['gradient'][label] = data.get('gradient', None)
                component_results_inv['hessian'][label] = data.get('hessian', None)

        # Actually analyze
        all_results = {}
        if all(v is not None for v in component_results_inv['energy'].values()):
            all_results['energy'] = self._analyze("energy", component_results_inv['energy'])
        if all(v is not None for v in component_results_inv['gradient'].values()):
            all_results['gradient'] = self._analyze(DriverEnum.gradient, component_results_inv['gradient'])
        if all(v is not None for v in component_results_inv['hessian'].values()):
            all_results['hessian'] = self._analyze(DriverEnum.hessian, component_results_inv['hessian'])

        return all_results