import abc
import copy

# printing and logging formatting niceties
import pprint
from functools import partial
import numpy as np
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
nppp = partial(np.array_str, max_line_width=120, precision=8, suppress_small=True)
nppp10 = partial(np.array_str, max_line_width=120, precision=10, suppress_small=True)

from ast import literal_eval
# v2: from typing import TYPE_CHECKING, Any, ClassVar, Dict, List, Tuple, Union, Literal, Optional
from typing import Any, Dict, List, Mapping, Tuple, Union, Literal, Optional

# v2: from pydantic import ConfigDict, field_validator, FieldValidationInfo, computed_field, BaseModel, Field
try:
    from pydantic.v1 import create_model, Field, validator, BaseModel
except ImportError:
    from pydantic import create_model, Field, validator, BaseModel

from qcelemental.models import FailedOperation, Molecule, DriverEnum, ProtoModel, AtomicResult, AtomicInput
import qcengine as qcng
from .manybody_v1 import BsseEnum, ManyBodyKeywords, ManyBodyInput, ManyBodyResult, ManyBodyResultProperties
from qcmanybody import ManyBodyCalculator
from qcmanybody.utils import delabeler, provenance_stamp


class BaseComputerQCNG(ProtoModel):
    """Base class for "computers" that plan, run, and process QC tasks."""

    @abc.abstractmethod
    def compute(self):
        pass

    @abc.abstractmethod
    def plan(self):
        pass

    # TODO can remove?
    # v2: model_config = ConfigDict(
    # v2:     extra="allow",
    # v2:     frozen=False,
    # v2: )
    class Config:
        extra = "allow"
        frozen = False


class AtomicComputerQCNG(BaseComputerQCNG):
    """Computer for analytic single-geometry computations."""

    molecule: Molecule = Field(..., description="The molecule to use in the computation.")
    basis: str = Field(..., description="The quantum chemistry basis set to evaluate (e.g., 6-31g, cc-pVDZ, ...).")
    method: str = Field(..., description="The quantum chemistry method to evaluate (e.g., B3LYP, MP2, ...).")
    driver: DriverEnum = Field(..., description="The resulting type of computation: energy, gradient, hessian, properties."
        "Note for finite difference that this should be the target driver, not the means driver.")
    keywords: Dict[str, Any] = Field(default_factory=dict, description="The keywords to use in the computation.")
    program: str = Field(..., description="Which program harness to run single-point with.")
    computed: bool = Field(False, description="Whether quantum chemistry has been run on this task.")
    result: Any = Field(default_factory=dict, description=":py:class:`~qcelemental.models.AtomicResult` return.")
    result_id: Optional[str] = Field(None, description="The optional ID for the computation.")

    # v2: @field_validator("basis")
    @validator("basis")
    @classmethod
    def set_basis(cls, basis):
        return basis.lower()

    # v2: @field_validator("method")
    @validator("method")
    @classmethod
    def set_method(cls, method):
        return method.lower()

    # v2: @field_validator("keywords")
    @validator("keywords")
    @classmethod
    def set_keywords(cls, keywords):
        return copy.deepcopy(keywords)

    def plan(self) -> AtomicInput:
        """Form QCSchema input from member data."""

        atomic_model = AtomicInput(**{
            "molecule": self.molecule,
            "driver": self.driver,
            "model": {
                "method": self.method,
                "basis": self.basis
            },
            "keywords": self.keywords,
        })

        return atomic_model

    def compute(self, client: Optional["qcportal.FractalClient"] = None) -> None:
        """Run quantum chemistry single-point.

        NOTE: client logic removed (compared to psi4.driver.AtomicComputer)
        """
        from qcengine import compute as qcng_compute

        if self.computed:
            return

        #logger.info(f'<<< JSON launch ... {self.molecule.name} {self.molecule.nuclear_repulsion_energy()}')

        self.result = qcng_compute(
            self.plan(),
            self.program,
            raise_error=False,  #True,
            #task_config=task_config,
        )

        #pp.pprint(self.result.model_dump())
        #logger.debug(pp.pformat(self.result.model_dump()))
        self.computed = True

    def get_results(self, client: Optional["qcportal.FractalClient"] = None) -> AtomicResult:
        """Return results as Atomic-flavored QCSchema.

        NOTE: client removed (compared to psi4.driver.AtomicComputer)
        """
        if self.result:
            return self.result


class ManyBodyComputerQCNG(BaseComputerQCNG):

    input_data: ManyBodyInput = Field(
        ...,
        description="Input schema containing the relevant settings for performing the many body "
            "expansion. This is entirely redundant with the piecemeal assembly of this Computer class "
            "and is only stored to be available for error handling and exact reconstruction of ManyBodyResult.",
    )
    bsse_type: List[BsseEnum] = Field(
        [BsseEnum.cp],
        # v2: description=ManyBodyKeywords.model_fields["bsse_type"].description,
        description=ManyBodyKeywords.__fields__["bsse_type"].field_info.description,
    )
    molecule: Molecule = Field(
        ...,
        description="Target molecule for many body expansion (MBE) or interaction energy (IE) analysis. "
            "Fragmentation should already be defined in `fragments` and related fields.",
    )
    driver: DriverEnum = Field(
        ...,
        description="The computation driver; i.e., energy, gradient, hessian. In case of ambiguity (e.g., MBE gradient "
            "through finite difference energies or MBE energy through composite method), this field refers to the "
            "*target* derivative, not any *means* specification.",
    )
    embedding_charges: Dict[int, List[float]] = Field(
        {},
        description="Atom-centered point charges to be used on molecule fragments whose basis sets are not included in "
            "the computation. Keys: 1-based index of fragment. Values: list of atom charges for that fragment.",
        json_schema_extra={
            "shape": ["nfr", "<varies: nat in ifr>"],
        },
    )
    return_total_data: Optional[bool] = Field(  # after driver, embedding_charges
        None,
        validate_default=True,
        # v2: description=ManyBodyKeywords.model_fields["return_total_data"].description,
        description=ManyBodyKeywords.__fields__["return_total_data"].field_info.description,
    )
    levels: Optional[Dict[Union[int, Literal["supersystem"]], str]] = Field(
        None,
        validate_default=True,
        # v2: description=ManyBodyKeywords.model_fields["levels"].description + \
        description=ManyBodyKeywords.__fields__["levels"].field_info.description + \
            "Examples above are processed in the ManyBodyComputer, and once processed, only the values should be used. "
            "The keys turn into nbodies_per_mc_level, as notated below. "
            "* {1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'} -> nbodies_per_mc_level=[[1], [2], ['supersystem']] "
            "* {2: 'ccsd(t)/cc-pvdz', 3: 'mp2'} -> nbodies_per_mc_level=[[1, 2], [3]] ",
    )
    max_nbody: Optional[int] = Field(
        None,
        validate_default=True,
        # v2: description=ManyBodyKeywords.model_fields["max_nbody"].description,
        description=ManyBodyKeywords.__fields__["max_nbody"].field_info.description,
    )
    supersystem_ie_only: Optional[bool] = Field(  # after max_nbody
        False,
        validate_default=True,
        # v2: description=ManyBodyKeywords.model_fields["supersystem_ie_only"].description,
        description=ManyBodyKeywords.__fields__["supersystem_ie_only"].field_info.description,
    )
    task_list: Dict[str, Any] = {}  #MBETaskComputers] = {}

    # TODO @computed_field(description="Number of distinct fragments comprising full molecular supersystem.")
    @property
    def nfragments(self) -> int:
        return len(self.molecule.fragments)

    # v2: @field_validator("bsse_type", mode="before")
    @validator("bsse_type", pre=True)
    @classmethod
    def set_bsse_type(cls, v: Any) -> List[BsseEnum]:
        if not isinstance(v, list):
            v = [v]
        # emulate ordered set
        # * bt.lower() as return (w/i `list(dict.fromkeys([bt.lower() ...`)
        #   works until aliases added to BsseEnum
        # * BsseEnum[bt].value as return works for good vals, but passing bad
        #   vals through as bt lets pydantic raise a clearer error message
        return list(dict.fromkeys([(BsseEnum[bt.lower()].value if bt.lower() in BsseEnum.__members__ else bt.lower()) for bt in v]))

    # v2: @field_validator("embedding_charges")
    @validator("embedding_charges", pre=True)
    @classmethod
    # v2: def set_embedding_charges(cls, v: Any, info: FieldValidationInfo) -> Dict[int, List[float]]:
    def set_embedding_charges(cls, v, values): # -> Dict[int, List[float]]:
        # v2: if len(v) != info.data["nfragments"]:
        if len(v) != values["nfragments"]:
            raise ValueError("embedding_charges dict should have entries for each 1-indexed fragment.")

        return v

    # v2: @field_validator("return_total_data")
    @validator("return_total_data", always=True)
    @classmethod
    # v2: def set_return_total_data(cls, v: Optional[bool], info: FieldValidationInfo) -> bool:
    def set_return_total_data(cls, v: Optional[bool], values) -> bool:
        print(f"hit return_total_data validator with {v}", end="")
        if v is not None:
            rtd = v
        # v2: elif info.data["driver"] in ["gradient", "hessian"]:
        elif values["driver"] in ["gradient", "hessian"]:
            rtd = True
        else:
            rtd = False

        # v2: if info.data.get("embedding_charges", False) and rtd is False:
        if values.get("embedding_charges", False) and rtd is False:
            raise ValueError("Cannot return interaction data when using embedding scheme.")

        print(f" ... setting rtd={rtd}")
        return rtd

    # v2: @field_validator("levels")
    @validator("levels", always=True)
    @classmethod
    # v2: def set_levels(cls, v: Any, info: FieldValidationInfo) -> Dict[Union[int, Literal["supersystem"]], str]:
    def set_levels(cls, v: Any, values) -> Dict[Union[int, Literal["supersystem"]], str]:
        print(f"hit levels validator with {v}", end="")

        if v is None:
            pass
            # TODO levels = {plan.max_nbody: method}
            #v = {info.data["nfragments"]: "???method???"}
            #v = {len(info.data["molecule"].fragments): "???method???"}
            # v2: v = {len(info.data["molecule"].fragments): "(auto)"}
            v = {len(values["molecule"].fragments): "(auto)"}
        else:
            # rearrange bodies in order with supersystem last lest body count fail in organization loop below
            v = dict(sorted(v.items(), key=lambda item: 1000 if item[0] == "supersystem" else item[0]))

            # rm 1 for cp-only
            # We define cp as being a correction to only interaction energies
            # If only doing cp, we need to ignore any user-specified 1st (monomer) level
            #if 'cp' in kwargs.get("bsse_type", None) and 'nocp' not in kwargs.get("bsse_type", None):
            #    if 1 in levels.keys():
            #        removed_level = levels.pop(1)
            #        logger.info("NOTE: User specified exclusively 'cp' correction, but provided level 1 details")
            #        logger.info(f"NOTE: Removing level {removed_level}")
            #        logger.info("NOTE: For total energies, add 'nocp' to bsse_list")

        print(f" ... setting levels={v}")
        return v

    # TODO @computed_field(
    # TODO     description="Distribution of active n-body levels among model chemistry levels. All bodies in range "
    # TODO         "[1, self.max_nbody] must be present exactly once. Number of items in outer list is how many different "
    # TODO         "modelchems. Each inner list specifies what n-bodies to be run at the corresponding modelchem (e.g., "
    # TODO         "`[[1, 2]]` has max_nbody=2 and 1-body and 2-body contributions computed at the same level of theory; "
    # TODO         "`[[1], [2]]` has max_nbody=2 and 1-body and 2-body contributions computed at different levels of theory. "
    # TODO         "An entry 'supersystem' means all higher order n-body effects up to the number of fragments. The n-body "
    # TODO         "levels are effectively sorted in the outer list, and any 'supersystem' element is at the end.")
        #json_schema_extra={
        #    "shape": ["nmc", "<varies>"],
        #},
    @property
    def nbodies_per_mc_level(self) -> List[List[Union[int, Literal["supersystem"]]]]:
        print(f"hit nbodies_per_mc_level", end="")

        # Organize nbody calculations into modelchem levels
        # * expand keys of `levels` into full lists of nbodies covered. save to plan, resetting max_nbody accordingly
        # * below, process values of `levels`, which are modelchem strings, into kwargs specs
        nbodies_per_mc_level = []
        prev_body = 0
        print("\nAAA", self.levels)
        for nb in self.levels:
            nbodies = []
            print("BBB bfore", nb, nbodies, prev_body)
            if nb == "supersystem":
                nbodies.append(nb)
            elif nb != (prev_body + 1):
                for m in range(prev_body + 1, nb + 1):
                    nbodies.append(m)
            else:
                nbodies.append(nb)
            print("BBB after", nb, nbodies)
            nbodies_per_mc_level.append(nbodies)
            prev_body = nb  # formerly buggy `+= 1`

        print(f" ... setting nbodies_per_mc_level={nbodies_per_mc_level}")
        return nbodies_per_mc_level

    # v2: @field_validator("max_nbody")
    @validator("max_nbody", always=True)
    @classmethod
    # v2: def set_max_nbody(cls, v: Any, info: FieldValidationInfo) -> int:
    def set_max_nbody(cls, v: Any, values) -> int:
        print(f"hit max_nbody validator with {v}", end="")
        # v2: levels_max_nbody = max(nb for nb in info.data["levels"] if nb != "supersystem")
        # v2: nfr = len(info.data["molecule"].fragments)
        levels_max_nbody = max(nb for nb in values["levels"] if nb != "supersystem")
        nfr = len(values["molecule"].fragments)
        print(f" {levels_max_nbody=} {nfr=}", end="")

        #ALT if v == -1:
        if v is None:
            v = levels_max_nbody
        elif v < 0 or v > nfr:
            raise ValueError(f"max_nbody={v} should be between 1 and {nfr}.")
        elif v != levels_max_nbody:
            #raise ValueError(f"levels={levels_max_nbody} contradicts user max_nbody={v}.")
            # TODO reconsider logic. move this from levels to here?
            # v2: info.data["levels"] = {v: "(auto)"}
            values["levels"] = {v: "(auto)"}
        else:
            pass
            # TODO once was           return min(v, nfragments)

        print(f" ... setting max_nbody={v}")
        return v

#       levels          max_nbody           F-levels        F-max_nbody     result
#
#       {stuff}         None                {stuff}         set from stuff  all consistent; max_nbody from levels
#       None            int                 {int: mtd}      int             all consistent; levels from max_nbody
#       None            None                {nfr: mtd}      nfr             all consistent; any order
#       {stuff}         int                 {stuff}         int             need to check consistency

    # TODO also, perhaps change nbodies_per_mc_level into dict of lists so that pos'n/label indexing coincides

    # v2: @field_validator("supersystem_ie_only")
    @validator("supersystem_ie_only", always=True)
    @classmethod
    # v2: def set_supersystem_ie_only(cls, v: Optional[bool], info: FieldValidationInfo) -> bool:
    def set_supersystem_ie_only(cls, v: Optional[bool], values) -> bool:
        print(f"hit supersystem_ie_only validator with {v}", end="")
        sio = v
        # v2: _nfr = len(info.data["molecule"].fragments)
        _nfr = len(values["molecule"].fragments)

        # v2: _max_nbody = info.data["max_nbody"]
        # get(..., None) b/c in v1, all fields processed even if max_nbody previously failed
        _max_nbody = values.get("max_nbody", None)
        if (sio is True) and (_max_nbody != _nfr):
            raise ValueError(f"Cannot skip intermediate n-body jobs when max_nbody={_max_nbody} != nfragments={_nfr}.")

        if (sio is True) and ("vmfc" in values["bsse_type"]):
            raise ValueError(f"Cannot skip intermediate n-body jobs when VMFC in bsse_type={values['bsse_type']}. Use CP instead.")

        print(f" ... setting {sio=}")
        return sio

    @classmethod
    def from_qcschema_ben(cls, input_model: ManyBodyInput):

        computer_model = cls(
            molecule=input_model.molecule,
            driver=input_model.specification.driver,
            # v2: **input_model.specification.keywords.model_dump(exclude_unset=True),
            **input_model.specification.keywords.dict(exclude_unset=True),
            input_data=input_model,  # storage, to reconstitute ManyBodyResult
        )
        nb_per_mc = computer_model.nbodies_per_mc_level

        comp_levels = {}
        for mc_level_idx, mtd in enumerate(computer_model.levels.values()):
            for lvl1 in nb_per_mc[mc_level_idx]:
                comp_levels[int(lvl1)] = mtd

        specifications = {}
        for mtd, spec in computer_model.input_data.specification.specification.items():
            spec = spec.dict()
            specifications[mtd] = {}
            specifications[mtd]["program"] = spec.pop("program")
            specifications[mtd]["specification"] = spec
            specifications[mtd]["specification"].pop("schema_name", None)
            specifications[mtd]["specification"].pop("protocols", None)

        calculator_cls = ManyBodyCalculator(
            computer_model.molecule,
            computer_model.bsse_type,
            comp_levels,
            computer_model.return_total_data,
            computer_model.supersystem_ie_only,
        )

        component_results = {}

        for chem, label, imol in calculator_cls.iterate_molecules():
            inp = AtomicInput(molecule=imol, **specifications[chem]["specification"])

            _, real, bas = delabeler(label)
            result = qcng.compute(inp, specifications[chem]["program"])

            if not result.success:
                print(result.error.error_message)
                raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

            # pull out stuff
            props = {"energy", "gradient", "hessian"}

            component_results[label] = {}

            for p in props:
                if hasattr(result.properties, f"return_{p}"):
                    v = getattr(result.properties, f"return_{p}")
                    if v is not None:
                        component_results[label][p] = v

        analyze_back = calculator_cls.analyze(component_results)
        analyze_back["nbody_number"] = len(component_results)

        return computer_model.get_results(external_results=analyze_back)

    @classmethod
    def from_qcschema(cls, input_model: ManyBodyInput, build_tasks: bool = False):

        computer_model = cls(
            molecule=input_model.molecule,
            driver=input_model.specification.driver,
            # v2: **input_model.specification.keywords.model_dump(exclude_unset=True),
            **input_model.specification.keywords.dict(exclude_unset=True),
            input_data=input_model,  # storage, to reconstitute ManyBodyResult
        )

        print("\n<<<  (Z) QCEngine harness ManyBodyComputerQCNG.from_qcschema  >>>")
        # v2: pprint.pprint(computer_model.model_dump(), width=200)
        pprint.pprint(computer_model.dict(), width=200)

        if build_tasks:
            # Note: self.input_data or input input_model.
            #   Also, we've always tried to keep build_tasks separate (see psi4.driver.task_planner).
            #   Clearly one _can_ init from ManyBodyInput -- make this sep fn?

            for mc_level_idx, mtd in enumerate(computer_model.levels.values()):
                atspec = computer_model.input_data.specification.specification[mtd]

                computer_model.build_tasks(
                    AtomicComputerQCNG,
                    mc_level_idx=mc_level_idx,
                    program=atspec.program,
                    driver=atspec.driver,  # use manybodyspecification.driver instead?
                    method=atspec.model.method,
                    basis=atspec.model.basis,
                    keywords=atspec.keywords,
                    #protocols=atspec.protocols,
                    #**kwargs,
                )

        return computer_model

    def build_tasks(
        self,
        mb_computer: AtomicComputerQCNG, #MBETaskComputers,
        mc_level_idx: int,
        **kwargs: Dict[str, Any],
    ) -> int:
        """Adds to the task_list as many new unique tasks as necessary to treat a single model chemistry level at one
        or several n-body levels. New tasks are of type *mb_computer* with model chemistry level specified in *kwargs*
        and n-body levels accessed through *mc_level_idx*.

        Parameters
        ----------
#        mb_computer
#            Class of task computers to instantiate and add to self.task_list. Usually :class:`~psi4.driver.AtomicComputer` but may be other when wrappers are layered.
#        mc_level_idx
#            Position in field self.nbodies_per_mc_level used to obtain ``nbodies``, the list of n-body
#            levels (e.g., `[1]` or `[1, 2]` or `["supersystem"]`) to which the modelchem specified in **kwargs** applies.
#            That is, `nbodies = self.nbodies_per_mc_level[mc_level_idx]`.
#            Note the natural 1-indexing of ``nbodies`` _contents_, so `[1]` covers one-body contributions.
#            The corresponding user label is the 1-indexed counterpart, `mc_level_lbl = mc_level_idx + 1`
#            Formerly nlevel as in `nbody = self.nbody_list[nbody_level=nlevel]`.
#        kwargs
#            Other arguments for initializing **mb_computer**. In particular, specifies model chemistry.

        Returns
        -------
        count : int
            Number of new tasks planned by this call.
            Formerly, didn't include supersystem in count.

        """
# qcng:        from psi4.driver.driver_nbody import build_nbody_compute_list

        # TODO method not coming from levels right

        # Get the n-body orders for this level. e.g., [1] or [2, 3] or ["supersystem"]
        nbodies = self.nbodies_per_mc_level[mc_level_idx]

#        info = "\n" + p4util.banner(f" ManyBody Setup: N-Body Levels {nbodies}", strNotOutfile=True) + "\n"
#        core.print_out(info)
#        logger.info(info)

#        for kwg in ['dft_functional']:
#            if kwg in kwargs:
#                kwargs['keywords']['function_kwargs'][kwg] = kwargs.pop(kwg)

        count = 0
        template = copy.deepcopy(kwargs)

        # Get compute list
        if nbodies == ["supersystem"]:
            # Add supersystem computation if requested -- always nocp
            data = template
            data["molecule"] = self.molecule
            key = f"supersystem_{self.nfragments}"
            self.task_list[key] = mb_computer(**data)
            count += 1

            compute_dict = build_nbody_compute_list(
                ["nocp"], list(range(1, self.max_nbody + 1)),
                self.nfragments, self.return_total_data, self.supersystem_ie_only)
        else:
            compute_dict = build_nbody_compute_list(
                self.bsse_type, nbodies,
                self.nfragments, self.return_total_data, self.supersystem_ie_only)

        def lab_labeler(item) -> str:
            # note 0-index to 1-index shift for label
            return f"{mc_level_idx + 1}_{item}"

        print("HHHH", compute_dict)
        # Add current compute list to the master task list
        # * `pair` looks like `((1,), (1, 3))` where first is real (not ghost) fragment indices
        #    and second is basis set fragment indices, all 1-indexed
        for nb in compute_dict["all"]:
            for pair in compute_dict["all"][nb]:
                lbl = lab_labeler(pair)
                if lbl in self.task_list:
                    continue

                data = template
                ghost = list(set(pair[1]) - set(pair[0]))
                # while psi4.core.Molecule.extract_subsets takes 1-indexed, qcel.models.Molecule.get_fragment takes 0-indexed.
                real0 = [(idx - 1) for idx in list(pair[0])]
                ghost0 = [(idx - 1) for idx in ghost]
                data["molecule"] = self.molecule.get_fragment(real=real0, ghost=ghost0, group_fragments=False)  # orient?
#                if self.embedding_charges:
#                    embedding_frags = list(set(range(1, self.nfragments + 1)) - set(pair[1]))
#                    charges = []
#                    for frag in embedding_frags:
#                        positions = self.molecule.extract_subsets(frag).geometry().np.tolist()
#                        charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[frag])])
#                    data['keywords']['function_kwargs'].update({'external_potentials': charges})

                self.task_list[lbl] = mb_computer(**data)
                count += 1

        return count

    def plan(self):
        # uncalled function
        return [t.plan() for t in self.task_list.values()]

    def compute(self, client: Optional["qcportal.FractalClient"] = None) -> None:
        """Run quantum chemistry.

        NOTE: client logic removed (compared to psi4.driver.ManyBodyComputer)
        """
# qcng:        from psi4.driver.p4util import banner

# qcng:        info = "\n" + banner(f" ManyBody Computations ", strNotOutfile=True) + "\n"
        #logger.info(info)

        for t in self.task_list.values():
            t.compute(client=client)

    def prepare_results(
        self,
        results: Optional[Dict[str, "MBETaskComputers"]] = None,
        client: Optional["qcportal.FractalClient"] = None,
    ) -> Dict[str, Any]:
        """Process the results from all n-body component molecular systems and model chemistry levels into final quantities.

        NOTE: client removed (compared to psi4.driver.ManyBodyComputer) (multilevel call, too)

#        Parameters
#        ----------
#        results
#            A set of tasks to process instead of self.task_list. Used in multilevel processing to pass a subset of
#            self.task_list filtered to only one modelchem level.
#        client
#            QCFractal client if using QCArchive for distributed compute.
#
#        Returns
#        -------
#        nbody_results
#            When the ManyBodyComputer specifies a single model chemistry level (see self.nbodies_per_mc_level), the
#            return is a dictionary, nbody_results, described in the table below. Many of the items are actually filled
#            by successive calls to assemble_nbody_components(). When multiple model chemistry levels are specified, this
#            function diverts its return to driver_nbody_multilevel.prepare_results() wherein each mc level calls this
#            function again and collects separate nbody_results dictionaries and processes them into a final return that
#            is a small subset of the table below.
#
#
#                                       ptype_size = (1,)/(nat, 3)/(3 * nat, 3 * nat)
#                                        e/g/h := energy or gradient or Hessian
#                                        rtd := return_total_data
#
        """
# qcng:        from psi4.driver.driver_nbody import assemble_nbody_components
# qcng:        from psi4.driver.driver_nbody_multilevel import prepare_results as multilevel_prepare_results

        if results is None:
            print(f"RESULTS setting empty")
            results = {}

#        # formerly nlevels
        mc_level_labels = {i.split("_")[0] for i in self.task_list}
        print("RESULTS", mc_level_labels, len(results), results.keys())
        if len(mc_level_labels) > 1 and not results:  # seeming fix to recursion
#        if len(self.nbodies_per_mc_level) > 1 and not results:  # max recursion
            print(f"RESULTS calling MLVL")
            return multilevel_prepare_results(self, client=client)

        results_list = {k: v.get_results(client=client) for k, v in (results.items() or self.task_list.items())}

#        pp.pprint(results_list["1_((1, 2), (1, 2))"].model_dump())
        #if "1_((1, 2, 3, 4), (1, 2, 3, 4))" in results_list:
        #    pp.pprint(results_list["1_((1, 2, 3, 4), (1, 2, 3, 4))"].model_dump())
        if "1_((3,), (3,))" in results_list:
            # v2: pp.pprint(results_list["1_((3,), (3,))"].model_dump())
            pp.pprint(results_list["1_((3,), (3,))"].dict())
        trove = {  # AtomicResult.properties return None if missing
            "energy": {k: v.properties.return_energy for k, v in results_list.items()},
            "gradient": {k: v.properties.return_gradient for k, v in results_list.items()},
            "hessian": {k: v.properties.return_hessian for k, v in results_list.items()},
        }

#        # TODO: make assemble_nbody_components and driver_nbody_multilevel.prepare_results into class functions.
#        #   note that the former uses metadata as read-only (except for one solveable case) while the latter overwrites self (!).
        metadata = {
            "quiet": False, #self.quiet,
            "nbodies_per_mc_level": self.nbodies_per_mc_level,
            "bsse_type": self.bsse_type,
            "nfragments": self.nfragments,
            "return_total_data": self.return_total_data,
            "supersystem_ie_only": self.supersystem_ie_only,
            "molecule": self.molecule,
            "embedding_charges": bool(self.embedding_charges),
            "max_nbody": self.max_nbody,
        }
        if self.driver.name == "energy":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())

        elif self.driver.name == "gradient":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())
            nbody_results.update(assemble_nbody_components("gradient", trove["gradient"], metadata.copy()))

        elif self.driver.name == "hessian":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())
            nbody_results.update(assemble_nbody_components("gradient", trove["gradient"], metadata.copy()))
            nbody_results.update(assemble_nbody_components("hessian", trove["hessian"], metadata.copy()))

        # save some mc_(frag, bas) component results
        # * formerly, intermediates_energy was intermediates2
        # * formerly, intermediates_gradient was intermediates_ptype
        # * formerly, intermediates_hessian was intermediates_ptype

        nbody_results["intermediates"] = {}
        for idx, task in results_list.items():
            mc, frag, bas = lab_delabeler(idx)
            nbody_results["intermediates"][f"N-BODY ({frag})@({bas}) TOTAL ENERGY"] = task.properties.return_energy

        nbody_results["intermediates_energy"] = trove["energy"]

        if not all(x is None for x in trove["gradient"].values()):
            nbody_results["intermediates_gradient"] = trove["gradient"]

        if not all(x is None for x in trove["hessian"].values()):
            nbody_results["intermediates_hessian"] = trove["hessian"]

#        debug = False
#        if debug:
#            for k, v in nbody_results.items():
#                if isinstance(v, np.ndarray):
#                    print(f"CLS-prepared results >>> {k} {v.size}")
#                elif isinstance(v, dict):
#                    print(f"CLS-prepared results >>> {k} {len(v)}")
#                    for k2, v2 in v.items():
#                        if isinstance(v2, np.ndarray):
#                            print(f"CLS-prepared results      >>> {k2} {v2.size}")
#                        else:
#                            print(f"CLS-prepared results      >>> {k2} {v2}")
#                else:
#                    print(f"CLS-prepared results >>> {k} {v}")

        return nbody_results


    def get_results(self, client: Optional["qcportal.FractalClient"] = None, external_results: Optional[Dict] = None) -> ManyBodyResult:
        """Return results as ManyBody-flavored QCSchema.

        NOTE: client removed (compared to psi4.driver.ManyBodyComputer)
        """
# qcng:        from psi4.driver.p4util import banner

# qcng:        info = "\n" + banner(f" ManyBody Results ", strNotOutfile=True) + "\n"
        #logger.info(info)

        if external_results is None:
            results = self.prepare_results(client=client)
            ret_energy = results.pop("ret_energy")
            ret_ptype = results.pop("ret_ptype")
            ret_gradient = results.pop("ret_gradient", None)
            nbody_number = len(self.task_list)
        else:
            ret_energy = external_results.pop("ret_energy")
            ret_ptype = ret_energy if self.driver == "energy" else external_results.pop(f"ret_{self.driver.name}")
            ret_gradient = external_results.pop("ret_gradient", None)
            nbody_number = external_results.pop("nbody_number")

        # load QCVariables
        qcvars = {
            'NUCLEAR REPULSION ENERGY': self.molecule.nuclear_repulsion_energy(),
            'NBODY NUMBER': nbody_number,
        }

        properties = {
            "calcinfo_nmc": len(self.nbodies_per_mc_level),
            "calcinfo_nfr": self.nfragments,  # or len(self.molecule.fragments)
            "calcinfo_natom": len(self.molecule.symbols),
            "calcinfo_nmbe": nbody_number,
            "nuclear_repulsion_energy": self.molecule.nuclear_repulsion_energy(),
            "return_energy": ret_energy,
        }

        if external_results is None:
            for k, val in results.items():
                qcvars[k] = val
        else:
            for k, val in external_results.items():
                if k == "results":
                    k = "nbody"
                qcvars[k] = val

        qcvars['CURRENT ENERGY'] = ret_energy
        if self.driver == 'gradient':
            qcvars['CURRENT GRADIENT'] = ret_ptype
            properties["return_gradient"] = ret_ptype
        elif self.driver == 'hessian':
            qcvars['CURRENT GRADIENT'] = ret_gradient
            qcvars['CURRENT HESSIAN'] = ret_ptype
            properties["return_gradient"] = ret_gradient
            properties["return_hessian"] = ret_ptype

#        build_out(qcvars)
        atprop = build_manybodyproperties(qcvars["nbody"])
        print("ATPROP")
        # v2: pp.pprint(atprop.model_dump())
        pp.pprint(atprop.dict())

#        output_data = {
#            "schema_version": 1,
#            "molecule": gamessmol,  # overwrites with outfile Cartesians in case fix_*=F
#            "extras": {**input_model.extras},
#            "native_files": {k: v for k, v in outfiles.items() if v is not None},
#            "properties": atprop,

#####
#        nbody_model = self.get_results(client=client)
#        ret = nbody_model.return_result
#
#        wfn = core.Wavefunction.build(self.molecule, "def2-svp", quiet=True)
#
#        # TODO all besides nbody may be better candidates for extras than qcvars. energy/gradient/hessian_body_dict in particular are too simple for qcvars (e.g., "2")

        print("QCVARS PRESCREEN")
        pp.pprint(qcvars)

        for qcv, val in qcvars.items():
            if isinstance(val, dict):
                if qcv in [
#            #"energies",  # retired
#            #"ptype",     # retired
#            "intermediates",
#            "intermediates_energy",  #"intermediates2",
#            "intermediates_gradient",  #"intermediates_ptype",
#            "intermediates_hessian",  #"intermediates_ptype",
#            "energy_body_dict",
#            "gradient_body_dict",  # ptype_body_dict
#            "hessian_body_dict",  # ptype_body_dict
#            "nbody",
#            "cp_energy_body_dict",
#            "nocp_energy_body_dict",
#            "vmfc_energy_body_dict",
#            "cp_gradient_body_dict",
#            "nocp_gradient_body_dict",
#            "vmfc_gradient_body_dict",
#            "cp_hessian_body_dict",
#            "nocp_hessian_body_dict",
#            "vmfc_hessian_body_dict",
                ]:
                    for qcv2, val2 in val.items():
#                        try:
                            qcvars[str(qcv2)] = val2
#                        except ValidationError:
#                            obj.set_variable(f"{self.driver.name} {qcv2}", val2)
            else:
                qcvars[qcv] = val

        # v2: component_results = self.model_dump()['task_list']  # TODO when/where include the indiv outputs
        component_results = self.dict()['task_list']  # TODO when/where include the indiv outputs
#        for k, val in component_results.items():
#            val['molecule'] = val['molecule'].to_schema(dtype=2)

        print("QCVARS")
        pp.pprint(qcvars)

        nbody_model = ManyBodyResult(
            **{
                'input_data': self.input_data,
                #'molecule': self.molecule,
                # v2: 'properties': {**atprop.model_dump(), **properties},
                'properties': {**atprop.dict(), **properties},
                'provenance': provenance_stamp(__name__),
                'extras': {
                    'qcvars': qcvars,
#                    'component_results': component_results,
                },
                'return_result': ret_ptype,
                'success': True,
            })

#        logger.debug('\nNBODY QCSchema:\n' + pp.pformat(nbody_model.model_dump()))

        return nbody_model





def lab_delabeler(item: str, return_obj: bool = False) -> Union[Tuple[str, str, str], Tuple[int, Tuple[int], Tuple[int]]]:
    """Transform labels like string "1_((2,), (1, 2))" into string tuple ("1", "2", "1, 2") or
    object tuple (1, (2,), (1, 2)).

    """
    mc, _, fragbas = item.partition("_")
    frag, bas = literal_eval(fragbas)

    if return_obj:
        return int(mc), frag, bas
    else:
        return mc, ", ".join(map(str, frag)), ", ".join(map(str, bas))


qcvars_to_manybodyproperties = {}
# v2: for skprop in ManyBodyResultProperties.model_fields.keys():
for skprop in ManyBodyResultProperties.__fields__.keys():
    qcvar = skprop.replace("_body", "-body").replace("_corr", "-corr").replace("_", " ").upper()
    qcvars_to_manybodyproperties[qcvar] = skprop
qcvars_to_manybodyproperties["CURRENT ENERGY"] = "return_energy"
qcvars_to_manybodyproperties["CURRENT GRADIENT"] = "return_gradient"
qcvars_to_manybodyproperties["CURRENT HESSIAN"] = "return_hessian"


def build_manybodyproperties(qcvars: Mapping) -> ManyBodyResultProperties:
    """For results extracted from QC output in QCDB terminology, translate to QCSchema terminology.

    Parameters
    ----------
    qcvars : PreservingDict
        Dictionary of calculation information in QCDB QCVariable terminology.

    Returns
    -------
    atprop : ManyBodyResultProperties
        Object of calculation information in QCSchema ManyBodyResultProperties terminology.

    """
    atprop = {}
    for pv, dpv in qcvars.items():
        if pv in qcvars_to_manybodyproperties:
            atprop[qcvars_to_manybodyproperties[pv]] = dpv

    return ManyBodyResultProperties(**atprop)
