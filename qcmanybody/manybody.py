from __future__ import annotations

import logging
import math
from typing import Set, Dict, Tuple, Union, Literal, Mapping, Any, Sequence

import numpy as np
from qcelemental.models import Molecule

from qcmanybody.assemble import sum_cluster_data
from qcmanybody.builder import build_nbody_compute_list
from qcmanybody.models import BsseEnum, FragBasIndex
from qcmanybody.utils import (
    delabeler,
    labeler,
    print_nbody_energy,
    collect_vars,
    copy_value,
    find_shape,
    shaped_zero,
    all_same_shape,
    expand_hessian,
    expand_gradient,
)

logger = logging.getLogger(__name__)


class ManyBodyCalculator:
    def __init__(
        self,
        molecule: Molecule,
        bsse_type: Sequence[BsseEnum],
        levels: Mapping[Union[int, Literal["supersystem"]], str],
        return_total_data: bool,
        supersystem_ie_only: bool,
    ):
        # TODO
        self.embedding_charges = {}

        self.molecule = molecule
        self.bsse_type = [BsseEnum(x) for x in bsse_type]
        self.return_total_data = return_total_data
        self.supersystem_ie_only = supersystem_ie_only
        self.nfragments = len(molecule.fragments)

        self.levels = levels

        # Levels without supersystem
        self.levels_no_ss = {int(k): v for k, v in levels.items() if k != "supersystem"}

        # Just a set of all the modelchems
        self.mc_levels = set(self.levels.values())

        self.max_nbody = max(self.levels_no_ss.keys())

        if len(self.bsse_type) == 0:
            raise ValueError("No BSSE correction specified")

        self.return_bsse_type = self.bsse_type[0]

        ###############################
        # Build nbodies_per_mc_level
        # TODO - use Lori's code
        # TODO - dict to list of lists to handle non-contiguous levels
        # TODO multilevel and supersystem_ie_only=T not allowed together
        # TODO supersystem in levels is not to be trusted -- nfrag only and skips levels
        max_level = max(self.levels_no_ss.keys())

        if set(range(1, max_level + 1)) != set(self.levels_no_ss.keys()):
            raise ValueError(f"Levels must be contiguous from 1 to {max_level}")

        self.nbodies_per_mc_level: Dict[str, list] = {mc_level: [] for mc_level in self.mc_levels}
        for k, v in self.levels_no_ss.items():
            self.nbodies_per_mc_level[v].append(k)

        self.nbodies_per_mc_level = {k: sorted(v) for k, v in self.nbodies_per_mc_level.items()}

        # Supersystem is always at the end
        if "supersystem" in levels:
            ss_mc = levels["supersystem"]
            self.nbodies_per_mc_level[ss_mc].append("supersystem")

        # To be built on the fly
        self.mc_compute_dict = None

        # Build size and slices dictionaries
        self.fragment_size_dict = {}
        self.fragment_slice_dict = {}
        iat = 0
        for ifr in range(1, self.nfragments + 1):
            nat = len(self.molecule.fragments[ifr - 1])
            self.fragment_size_dict[ifr] = nat
            self.fragment_slice_dict[ifr] = slice(iat, iat + nat)
            iat += nat

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
                self.supersystem_ie_only,
                self.max_nbody,
            )

        return self.mc_compute_dict

    def expand_gradient(self, grad: np.ndarray, bas: Tuple[int, ...]) -> np.ndarray:
        return expand_gradient(grad, bas, self.fragment_size_dict, self.fragment_slice_dict)

    def expand_hessian(self, hess: np.ndarray, bas: Tuple[int, ...]) -> np.ndarray:
        return expand_hessian(hess, bas, self.fragment_size_dict, self.fragment_slice_dict)

    def iterate_molecules(self) -> Tuple[str, str, Molecule]:
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
                    mol = self.molecule.get_fragment(real_atoms_0, ghost_atoms_0, orient=False, group_fragments=False)
                    mol = mol.copy(update={"fix_com": True, "fix_orientation": True})

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
        property_label: str,
        component_results: Dict[str, Union[float, np.ndarray]],
    ) -> Dict[str, Any]:
        """Assembles N-body components for a single derivative level and a single model chemistry level
        into interaction quantities according to requested BSSE treatment(s).
        """

        # which level are we assembling?
        delabeled = [delabeler(k) for k in component_results.keys()]
        mc_level_labels = {x[0] for x in delabeled}

        if len(mc_level_labels) != 1:
            raise RuntimeError(f"Multiple model chemistries passed into _assemble_nbody_components: {mc_level_labels}")

        mc_level = mc_level_labels.pop()
        if mc_level not in self.mc_levels:
            raise RuntimeError(f"Model chemistry {mc_level} not found in {self.mc_levels}")

        # get the range of nbodies and the required calculations for this level
        bsse_type = self.bsse_type
        return_bsse_type = self.return_bsse_type
        nbodies = self.nbodies_per_mc_level[mc_level]
        if "supersystem" in nbodies:
            nbodies = list(range(1, self.max_nbody + 1))
            bsse_type = [BsseEnum.nocp]
            return_bsse_type = BsseEnum.nocp

        max_nbody = max(nbodies)
        compute_dict = self.compute_map[mc_level]

        if not all_same_shape(component_results.values()):
            raise ValueError("All values in data dictionary must have the same shape.")

        # Use first data value to determine shape
        first_key = next(iter(component_results.keys()))
        property_shape = find_shape(component_results[first_key])

        # Final dictionaries
        cp_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        nocp_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        vmfc_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}

        cp_body_dict = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        nocp_body_dict = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        vmfc_body_dict = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}

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
            cp_by_level[nb] = sum_cluster_data(component_results, cp_compute_list[nb], mc_level)
            nocp_by_level[nb] = sum_cluster_data(component_results, nocp_compute_list[nb], mc_level)
            if nb in compute_dict["vmfc_levels"]:
                vmfc_by_level[nb] = sum_cluster_data(
                    component_results, compute_dict["vmfc_levels"][nb], mc_level, vmfc=True, nb=nb
                )

        # Extract data for monomers in monomer basis for CP total data
        if 1 in nbodies:
            monomers_in_monomer_basis = [v for v in compute_dict["nocp"][1] if len(v[1]) == 1]
            monomer_sum = sum_cluster_data(component_results, set(monomers_in_monomer_basis), mc_level)
        else:
            monomer_sum = shaped_zero(property_shape)

        # Compute cp
        if BsseEnum.cp in bsse_type:
            for nb in range(1, nbodies[-1] + 1):
                if nb == self.nfragments:
                    cp_body_dict[nb] = cp_by_level[nb] - bsse
                    continue

                for k in range(1, nb + 1):
                    take_nk = math.comb(self.nfragments - k - 1, nb - k)
                    sign = (-1) ** (nb - k)
                    cp_body_dict[nb] += take_nk * sign * cp_by_level[k]

                if nb == 1:
                    bsse = cp_body_dict[nb] - monomer_sum
                    cp_body_dict[nb] = copy_value(monomer_sum)
                else:
                    cp_body_dict[nb] -= bsse

        # Compute nocp
        if BsseEnum.nocp in bsse_type:
            for nb in range(1, nbodies[-1] + 1):
                if nb == self.nfragments:
                    nocp_body_dict[nb] = nocp_by_level[nb]
                    continue

                for k in range(1, nb + 1):
                    take_nk = math.comb(self.nfragments - k - 1, nb - k)
                    sign = (-1) ** (nb - k)
                    nocp_body_dict[nb] += take_nk * sign * nocp_by_level[k]

        # Compute vmfc
        if BsseEnum.vmfc in bsse_type:
            for nb in nbodies:
                # TODO I think this is correct for all properties...
                # for k in range(1, nb + 1):
                #    vmfc_body_dict[nb] += vmfc_by_level[k]

                # TODO - but below was used for gradient/hessian in psi4?
                if property_label == "energy":
                    for k in range(1, nb + 1):
                        vmfc_body_dict[nb] += vmfc_by_level[k]
                else:
                    if nb > 1:
                        vmfc_body_dict[nb] = vmfc_by_level[nb - 1]
                    vmfc_body_dict[nb] += vmfc_by_level[nb]

        # Collect specific and generalized returns
        results = {
            f"cp_{property_label}_body_dict": cp_body_dict,
            f"nocp_{property_label}_body_dict": nocp_body_dict,
            f"vmfc_{property_label}_body_dict": vmfc_body_dict,
        }

        # Overall return body dict & value for this property
        results[f"{property_label}_body_dict"] = results[f"{return_bsse_type.value}_{property_label}_body_dict"]
        results[f"ret_{property_label}"] = copy_value(results[f"{property_label}_body_dict"][max_nbody])

        if not self.return_total_data:
            results[f"ret_{property_label}"] -= results[f"{property_label}_body_dict"][1]

        return results

    def _analyze(
        self,
        property_label: str,
        property_results: Dict[str, Union[float, np.ndarray]],  # Label to results
    ):
        # Initialize with zeros
        if not all_same_shape(property_results.values()):
            raise ValueError("All values in data dictionary must have the same shape.")

        # Use first data value to determine shape
        first_key = next(iter(property_results.keys()))
        property_shape = find_shape(property_results[first_key])

        property_result = shaped_zero(property_shape)
        property_body_dict = {b.value: {} for b in self.bsse_type}
        property_body_contribution = {b.value: {} for b in self.bsse_type}

        # results per model chemistry
        mc_results = {}

        # sort by nbody level, ignore supersystem
        sorted_nbodies = [(k, v) for k, v in self.nbodies_per_mc_level.items() if v != ["supersystem"]]
        sorted_nbodies = sorted(sorted_nbodies, reverse=True, key=lambda x: x[1])
        for mc_label, nbody_list in sorted_nbodies:
            # filter to only one model chemistry
            filtered_results = {k: v for k, v in property_results.items() if delabeler(k)[0] == mc_label}

            if not filtered_results:
                if nbody_list == [1]:
                    # Note A.2: Note A.1 holds, but for the special case of CP-only
                    #   and rtd=False and multilevel with a separate level for
                    #   1-body, the skipped tasks run afoul of sanity checks, so
                    #   we'll add a dummy result.
                    filtered_results = {labeler(mc_label, [1000], [1000]): shaped_zero(property_shape)}
                else:
                    raise RuntimeError(f"No data found for model chemistry {mc_label}")

            nb_component_results = self._assemble_nbody_components(property_label, filtered_results)
            mc_results[mc_label] = nb_component_results

            for n in nbody_list[::-1]:
                property_bsse_dict = {b.value: shaped_zero(property_shape) for b in self.bsse_type}

                for m in range(n - 1, n + 1):
                    if m == 0:
                        continue

                    # Subtract the (n-1)-body contribution from the n-body contribution to get the n-body effect
                    sign = (-1) ** (1 - m // n)
                    for b in self.bsse_type:
                        property_bsse_dict[b.value] += (
                            sign * mc_results[mc_label][f"{b.value}_{property_label}_body_dict"][m]
                        )

                property_result += property_bsse_dict[self.return_bsse_type]
                for b in self.bsse_type:
                    property_body_contribution[b.value][n] = property_bsse_dict[b.value]

        if self.has_supersystem:
            # Get the MC label for supersystem tasks
            supersystem_mc_level = self.levels["supersystem"]

            # Super system recovers higher order effects at a lower level
            frag_range = tuple(range(1, self.nfragments + 1))

            ss_cresults = {k: v for k, v in property_results.items() if delabeler(k)[0] == supersystem_mc_level}
            ss_component_results = self._assemble_nbody_components(property_label, ss_cresults)
            mc_results[supersystem_mc_level] = ss_component_results

            # Compute components at supersystem level of theory
            ss_label = labeler(supersystem_mc_level, frag_range, frag_range)
            supersystem_result = property_results[ss_label]
            property_result += supersystem_result - ss_component_results[f"{property_label}_body_dict"][self.max_nbody]

            for b in self.bsse_type:
                property_body_contribution[b][self.nfragments] = (
                    supersystem_result - ss_component_results[f"{property_label}_body_dict"][self.max_nbody]
                )

        for b in self.bsse_type:
            bstr = b.value
            for n in property_body_contribution[bstr]:
                property_body_dict[bstr][n] = sum(
                    [
                        property_body_contribution[bstr][i]
                        for i in range(1, n + 1)
                        if i in property_body_contribution[bstr]
                    ]
                )

        if not self.return_total_data:
            # Remove monomer contribution for interaction data
            property_result -= property_body_dict[self.return_bsse_type][1]

        nbody_results = {
            f"ret_{property_label}": property_result,
            f"{property_label}_body_dict": property_body_dict,
            "mc_results": mc_results,
        }
        return nbody_results

    def analyze(
        self,
        component_results: Dict[str, Dict[str, Union[float, np.ndarray]]],  # component_results[label][property] = 1.23
    ):

        # All properties that were passed to us
        available_properties = set()
        for property_data in component_results.values():
            available_properties.update(property_data.keys())

        # reorganize to component_results_inv[property][label] = 1.23
        component_results_inv = {k: {} for k in available_properties}

        for cluster_label, property_data in component_results.items():
            for property_label, property_value in property_data.items():
                component_results_inv[property_label][cluster_label] = property_value

        # Remove any missing data
        component_results_inv = {k: v for k, v in component_results_inv.items() if v}

        # Actually analyze
        all_results = {}

        for property_label, property_results in component_results_inv.items():
            # Expand gradient and hessian
            if property_label == "gradient":
                property_results = {k: self.expand_gradient(v, delabeler(k)[2]) for k, v in property_results.items()}
            if property_label == "hessian":
                property_results = {k: self.expand_hessian(v, delabeler(k)[2]) for k, v in property_results.items()}

            r = self._analyze(property_label, property_results)
            all_results.update(r)

        # Analyze the total results
        nbody_dict = {}

        is_embedded = bool(self.embedding_charges)

        for b in self.bsse_type:
            print_nbody_energy(
                all_results["energy_body_dict"][b],
                f"{b.upper()}-corrected multilevel many-body expansion",
                self.nfragments,
                is_embedded,
            )

            if not self.has_supersystem:  # skipped levels?
                nbody_dict.update(
                    collect_vars(b.upper(), all_results["energy_body_dict"][b], self.max_nbody, is_embedded, self.supersystem_ie_only)
                )

        all_results["results"] = nbody_dict

        # Make dictionary with "1cp", "2cp", etc
        ebd = all_results["energy_body_dict"]
        all_results["energy_body_dict"] = {str(k) + b: v for b in ebd for k, v in ebd[b].items()}

        return all_results
