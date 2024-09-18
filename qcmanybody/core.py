from __future__ import annotations

import logging
import math
import os
import string
from collections import Counter, defaultdict
from typing import Any, Dict, Iterable, Literal, Mapping, Sequence, Set, Tuple, Union

import numpy as np
from qcelemental.models import Molecule

from qcmanybody.builder import build_nbody_compute_list
from qcmanybody.models import BsseEnum, FragBasIndex
from qcmanybody.utils import (
    all_same_shape,
    collect_vars,
    copy_value,
    delabeler,
    find_shape,
    labeler,
    modelchem_labels,
    print_nbody_energy,
    resize_gradient,
    resize_hessian,
    shaped_zero,
    sum_cluster_data,
)

logger = logging.getLogger(__name__)


__all__ = ["ManyBodyCalculator", "ManyBodyCore"]


class ManyBodyCore:
    def __init__(
        self,
        molecule: Molecule,
        bsse_type: Sequence[BsseEnum],
        levels: Mapping[Union[int, Literal["supersystem"]], str],
        *,
        return_total_data: bool,
        supersystem_ie_only: bool,
        embedding_charges: Mapping[int, Sequence[float]],
    ):
        self.embedding_charges = embedding_charges
        if self.embedding_charges:
            if not bool(os.environ.get("QCMANYBODY_EMBEDDING_CHARGES", False)):  # obscure until further validation
                raise ValueError(
                    f"Embedding charges for EE-MBE are still in testing. Set environment variable QCMANYBODY_EMBEDDING_CHARGES=1 to use at your own risk."
                )

        if isinstance(molecule, dict):
            mol = Molecule(**molecule)
        elif isinstance(molecule, Molecule):
            mol = molecule.copy()
        else:
            raise ValueError(f"Molecule input type of {type(molecule)} not understood.")
        self.molecule = mol
        self.bsse_type = [BsseEnum(x) for x in bsse_type]
        self.return_total_data = return_total_data
        self.supersystem_ie_only = supersystem_ie_only
        self.nfragments = len(self.molecule.fragments)

        self.levels = levels

        # Levels without supersystem
        self.levels_no_ss = {int(k): v for k, v in levels.items() if k != "supersystem"}

        # Just a set of all the modelchems
        self.mc_levels = set(self.levels.values())

        self.max_nbody = max(self.levels_no_ss.keys())

        if len(self.bsse_type) == 0:
            raise ValueError("No BSSE correction specified")

        if BsseEnum.vmfc in self.bsse_type and len(set(self.levels.values())) == 1:
            # For single-modelchem VMFC, NOCP & sometimes CP are produced for free
            if BsseEnum.nocp not in self.bsse_type:
                self.bsse_type.append(BsseEnum.nocp)
            if BsseEnum.cp not in self.bsse_type and self.max_nbody == self.nfragments:
                self.bsse_type.append(BsseEnum.cp)

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

        # order nbodies_per_mc_level keys (modelchems) by the lowest n-body level covered; any
        #   supersystem key (replaced below) is at the end. Order nbodies within each modelchem.
        #   Reset mc_levels to match.
        self.nbodies_per_mc_level = {
            k: sorted(v)
            for (k, v) in sorted(self.nbodies_per_mc_level.items(), key=lambda item: sorted(item[1] or [1000])[0])
        }
        assert self.mc_levels == set(self.nbodies_per_mc_level.keys())  # remove after some downstream testing
        self.mc_levels = self.nbodies_per_mc_level.keys()

        for mc, nbs in self.nbodies_per_mc_level.items():
            if nbs and ((nbs[-1] - nbs[0]) != len(nbs) - 1):
                raise ValueError(
                    f"QCManyBody: N-Body levels must be contiguous within a model chemistry spec ({mc}: {nbs}). Use an alternate spec name to accomplish this input."
                )
                # TODO - test and reenable if appropriate. my guess is that noncontig nb is fine on the core computing side,
                #   but trouble for computer and nbodies_per_mc_level inverting and indexing. Safer to deflect for now since input tweak allows the calc.

        # Supersystem is always at the end
        if "supersystem" in levels:
            ss_mc = levels["supersystem"]
            self.nbodies_per_mc_level[ss_mc].append("supersystem")

        # To be built on the fly
        self.mc_compute_dict = None

        if self.nfragments == 1:
            # Usually we try to "pass-through" edge cases, so a single-fragment mol would return 0 or ordinary energy,
            #   depending on rtd=T/F. But it seems more likely that user just forgot the fragments field, so we don't
            #   want to start a full energy on monsterMol. Reconsider handling in future.
            raise ValueError("""QCManyBody: Molecule fragmentation has not been specified through `fragments` field.""")

        if not np.array_equal(np.concatenate(self.molecule.fragments), np.arange(len(self.molecule.symbols))):
            raise ValueError("""QCManyBody: non-contiguous fragments could be implemented but aren't at present""")

        # Build size and slices dictionaries. Assumes fragments are contiguous
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

    def format_calc_plan(self, sset: str = "all") -> Tuple[str, Dict[str, Dict[int, int]]]:
        """Formulate per-modelchem and per-body job count data and summary text.

        Parameters
        ----------
        sset
            Among {"all", "nocp", "cp", "vmfc_compute"}, which data structure to return.

        Returns
        -------
        info
            A text summary with per- model chemistry and per- n-body-level job counts.
            ```
            Model chemistry "c4-ccsd" (§A):         22
                 Number of 1-body computations:     16 (nocp: 0, cp: 0, vmfc_compute: 16)
                 Number of 2-body computations:      6 (nocp: 0, cp: 0, vmfc_compute: 6)

            Model chemistry "c4-mp2" (§B):          28
                 Number of 1-body computations:     12 (nocp: 0, cp: 0, vmfc_compute: 12)
                 Number of 2-body computations:     12 (nocp: 0, cp: 0, vmfc_compute: 12)
                 Number of 3-body computations:      4 (nocp: 0, cp: 0, vmfc_compute: 4)
            ```
        Dict[str, Dict[int, int]]
            Data structure with outer key mc-label, inner key 1-indexed n-body, and value job count.
        """
        # Rearrange compute_list from key nb having values (species) to compute all of that nb
        #   to key nb having values counting that nb.
        compute_list_count = {}
        for mc, compute_dict in self.compute_map.items():
            compute_list_count[mc] = {}
            for sub in compute_dict:  # all, nocp, cp, vmfc
                all_calcs = set().union(*compute_dict[sub].values())
                compute_list_count[mc][sub] = Counter([len(frag) for (frag, _) in all_calcs])

        mc_labels = modelchem_labels(self.nbodies_per_mc_level, presorted=True)
        full_to_ordinal_mc_lbl = {v[0]: v[1] for v in mc_labels.values()}
        info = []
        for mc, counter in compute_list_count.items():
            all_counter = counter["all"]
            mcheader = f'    Model chemistry "{mc}" ({full_to_ordinal_mc_lbl[mc]}):'
            info.append(f"{mcheader:38} {sum(all_counter.values()):6}")
            for nb, count in sorted(all_counter.items()):
                other_counts = [f"{sub}: {counter[sub][nb]}" for sub in ["nocp", "cp", "vmfc_compute"]]
                info.append(f"        Number of {nb}-body computations: {count:6} ({', '.join(other_counts)})")
            info.append("")
        info = "\n".join(info)

        logger.info(info)
        return info, {mc: dsset[sset] for mc, dsset in compute_list_count.items()}

    def resize_gradient(self, grad: np.ndarray, bas: Tuple[int, ...], *, reverse: bool = False) -> np.ndarray:
        return resize_gradient(grad, bas, self.fragment_size_dict, self.fragment_slice_dict, reverse=reverse)

    def resize_hessian(self, hess: np.ndarray, bas: Tuple[int, ...], *, reverse: bool = False) -> np.ndarray:
        return resize_hessian(hess, bas, self.fragment_size_dict, self.fragment_slice_dict, reverse=reverse)

    def iterate_molecules(self) -> Iterable[Tuple[str, str, Molecule]]:
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

                    if self.embedding_charges:
                        embedding_frags = list(set(range(1, self.nfragments + 1)) - set(basis_atoms))
                        charges = []
                        for ifr in embedding_frags:
                            positions = self.molecule.get_fragment(ifr - 1).geometry.tolist()
                            charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[ifr])])
                        mol.extras["embedding_charges"] = charges

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

        # Accumulation dictionaries
        # * {bsse_type}_by_level is filled by sum_cluster_data to contain for NOCP
        #   & CP the summed total energies (or other property) of each nb-body. That is:
        #   * NOCP: {1: 1b@1b,    2: 2b@2b,      ..., max_nbody: max_nbody-b@max_nbody-b} and
        #   * CP:   {1: 1b@nfr-b, 2: 2b@nfr-b,   ..., max_nbody: max_nbody-b@nfr-b}.
        #   VMFC bookkeeping is different. For key 1 it contains the summed 1b total energies.
        #   But for higher keys, it contains each nb-body (non-additive) contribution to the energy.
        #   * VMFC: {1: 1b@1b,    2: 2b contrib, ..., max_nbody: max_nbody-b contrib}
        cp_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        nocp_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        vmfc_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}

        # * {bsse_type}_body_dict is usually filled with total energies (or other property).
        #   Multiple model chemistry levels may be involved.
        #   Generally, all consecutive keys between 1 and max_nbody will be present in the body_dict,
        #   but if supersystem_ie_only=T, only 1b and nfr-b are present, or if "supersystem" in levels, ???
        #   * TOT: {1: 1b@1b, 2: 2b tot prop with bsse_type treatment, ..., max_nbody: max_nbody-b tot prop with bsse_type treatment}
        #   If 1b@1b (monomers in monomer basis) aren't available, which can happen when return_total_data=F
        #   and 1b@1b aren't otherwise needed, body_dict contains interaction energies (or other property).
        #   * IE: {1: shaped_zero, 2: 2b interaction prop using bsse_type, ..., max_nbody: max_nbody-b interaction prop using bsse_type}
        #   For both TOT and IE cases, body_dict values are cummulative, not additive. For TOT, total,
        #   interaction, and contribution data in ManyBodyResultProperties can be computed in
        #   collect_vars. For IE, interaction and contribution data can be computed.
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
                for k in range(1, nb + 1):
                    vmfc_body_dict[nb] += vmfc_by_level[k]

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
        property_body_dict = {bt.value: {} for bt in self.bsse_type}
        property_body_contribution = {bt.value: {} for bt in self.bsse_type}

        # results per model chemistry
        mc_results = {}
        species_results = {}

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
                property_bsse_dict = {bt.value: shaped_zero(property_shape) for bt in self.bsse_type}

                for m in range(n - 1, n + 1):
                    if m == 0:
                        continue

                    # Subtract the (n-1)-body contribution from the n-body contribution to get the n-body effect
                    sign = (-1) ** (1 - m // n)
                    for bt in self.bsse_type:
                        property_bsse_dict[bt.value] += (
                            sign * mc_results[mc_label][f"{bt.value}_{property_label}_body_dict"][m]
                        )

                property_result += property_bsse_dict[self.return_bsse_type]
                for bt in self.bsse_type:
                    property_body_contribution[bt.value][n] = property_bsse_dict[bt.value]

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

            for bt in self.bsse_type:
                property_body_contribution[bt][self.nfragments] = (
                    supersystem_result - ss_component_results[f"{property_label}_body_dict"][self.max_nbody]
                )

        for bt in self.bsse_type:
            bstr = bt.value
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
        component_results: Dict[str, Dict[str, Union[float, np.ndarray]]],
    ):
        """

        Parameters
        ----------
        component_results
            Nested dictionary with results from all individual molecular system
            calculations, including all subsystem/basis combinations, all model
            chemistries, and all properties (e.g., e/g/h).

            For example, the below is the format for a nocp gradient run on a
            helium dimer with 1-body at CCSD and 2-body at MP2. The outer string
            key can be generated with the ``qcmanybody.utils.labeler`` function.
            The inner string key is any property; QCManyBody presently knows how
            to process energy/gradient/Hessian.
            ```
            {'["ccsd", [1], [1]]': {'energy': -2.87, 'gradient': array([[0., 0., 0.]])},
             '["ccsd", [2], [2]]': {'energy': -2.87, 'gradient': array([[0., 0., 0.]])},
             '["mp2", [1], [1]]': {'energy': -2.86, 'gradient': array([[0., 0., 0.]])},
             '["mp2", [2], [2]]': {'energy': -2.86, 'gradient': array([[0., 0., 0.]])},
             '["mp2", [1, 2], [1, 2]]': {'energy': -5.73, 'gradient': array([[ 0., 0., 0.0053], [ 0., 0., -0.0053]])},
            }
            ```
        """

        # All properties that were passed to us
        # * seed with "energy" so free/no-op jobs can process
        available_properties: Set[str] = {"energy"}
        for property_data in component_results.values():
            available_properties.update(property_data.keys())

        # reorganize to component_results_inv[property][label] = 1.23
        component_results_inv = {k: {} for k in available_properties}

        for cluster_label, property_data in component_results.items():
            for property_label, property_value in property_data.items():
                component_results_inv[property_label][cluster_label] = property_value

        # Remove any missing data
        component_results_inv = {k: v for k, v in component_results_inv.items() if v}
        if not component_results_inv:
            # Note B: Rarely, "no results" is expected, like for CP-only,
            #   rtd=False, and max_nbody=1. We'll add a dummy entry so
            #   processing can continue.
            component_results_inv["energy"] = {'["dummy", [1000], [1000]]': 0.0}

        # Actually analyze
        is_embedded = bool(self.embedding_charges)
        component_properties = defaultdict(dict)
        all_results = {}
        nbody_dict = {}
        stdout = ""
        #        all_results["energy_body_dict"] = {"cp": {1: 0.0}}

        for property_label, property_results in component_results_inv.items():
            # Expand gradient and hessian
            if property_label == "gradient":
                property_results = {k: self.resize_gradient(v, delabeler(k)[2]) for k, v in property_results.items()}
            if property_label == "hessian":
                property_results = {k: self.resize_hessian(v, delabeler(k)[2]) for k, v in property_results.items()}

            r = self._analyze(property_label, property_results)
            for k, v in property_results.items():
                component_properties[k]["calcinfo_natom"] = len(self.molecule.symbols)
                component_properties[k][f"return_{property_label}"] = v
            all_results.update(r)

        for bt in self.bsse_type:
            stdout += print_nbody_energy(
                all_results["energy_body_dict"][bt],
                f"{bt.formal()} ({bt.abbr()})",
                self.nfragments,
                modelchem_labels(self.nbodies_per_mc_level, presorted=True),
                is_embedded,
                self.supersystem_ie_only,
                self.max_nbody if self.has_supersystem else None,
            )

        for property_label in available_properties:
            for bt in self.bsse_type:
                nbody_dict.update(
                    collect_vars(
                        bt,
                        property_label,
                        all_results[f"{property_label}_body_dict"][bt],
                        self.max_nbody,
                        is_embedded,
                        self.supersystem_ie_only,
                        self.has_supersystem,
                    )
                )

        all_results["results"] = nbody_dict
        all_results["component_properties"] = component_properties

        # Make dictionary with "1cp", "2cp", etc
        ebd = all_results["energy_body_dict"]
        all_results["energy_body_dict"] = {str(k) + bt: v for bt in ebd for k, v in ebd[bt].items()}
        all_results["stdout"] = stdout

        return all_results


class ManyBodyCalculator(ManyBodyCore):
    # retire after a grace period
    def __init__(
        self,
        molecule: Molecule,
        bsse_type: Sequence[BsseEnum],
        levels: Mapping[Union[int, Literal["supersystem"]], str],
        return_total_data: bool,
        supersystem_ie_only: bool,
        embedding_charges: Mapping[int, Sequence[float]],
    ):
        super().__init__(
            molecule,
            bsse_type,
            levels,
            return_total_data=return_total_data,
            supersystem_ie_only=supersystem_ie_only,
            embedding_charges=embedding_charges,
        )
