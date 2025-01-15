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
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Mapping, Optional, Tuple, Union

# v2: from pydantic import ConfigDict, field_validator, FieldValidationInfo, computed_field, BaseModel, Field
try:
    from pydantic.v1 import BaseModel, Field, create_model, validator
except ImportError:
    from pydantic import create_model, Field, validator, BaseModel

from qcelemental.models import AtomicInput, AtomicResult, DriverEnum, FailedOperation, Molecule, ProtoModel

from qcmanybody import ManyBodyCore
from qcmanybody.models import BsseEnum, ManyBodyInput, ManyBodyKeywords, ManyBodyResult, ManyBodyResultProperties
from qcmanybody.utils import delabeler, provenance_stamp

if TYPE_CHECKING:
    import qcportal


__all__ = ["ManyBodyComputer"]


class BaseComputerQCNG(ProtoModel):
    """Base class for "computers" that plan, run, and process QC tasks."""

    @abc.abstractmethod
    def compute(self):
        pass

    @abc.abstractmethod
    def plan(self):
        pass

    # TODO can remove?
    #
    # v2: model_config = ConfigDict(
    # v2:     extra="allow",
    # v2:     frozen=False,
    # v2: )
    class Config:
        extra = "allow"
        # v2: frozen = False
        allow_mutation = True


class AtomicComputer(BaseComputerQCNG):
    """Computer for analytic single-geometry computations."""

    molecule: Molecule = Field(..., description="The molecule to use in the computation.")
    basis: str = Field(..., description="The quantum chemistry basis set to evaluate (e.g., 6-31g, cc-pVDZ, ...).")
    method: str = Field(..., description="The quantum chemistry method to evaluate (e.g., B3LYP, MP2, ...).")
    driver: DriverEnum = Field(
        ...,
        description="The resulting type of computation: energy, gradient, hessian, properties."
        "Note for finite difference that this should be the target driver, not the means driver.",
    )
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

        atomic_model = AtomicInput(
            **{
                "molecule": self.molecule,
                "driver": self.driver,
                "model": {"method": self.method, "basis": self.basis},
                "keywords": self.keywords,
            }
        )

        return atomic_model

    def compute(self, client: Optional["qcportal.FractalClient"] = None) -> None:
        """Run quantum chemistry single-point.

        NOTE: client logic removed (compared to psi4.driver.AtomicComputer)
        """
        from qcengine import compute as qcng_compute

        if self.computed:
            return

        # logger.info(f'<<< JSON launch ... {self.molecule.name} {self.molecule.nuclear_repulsion_energy()}')

        self.result = qcng_compute(
            self.plan(),
            self.program,
            raise_error=False,  # True,
            # task_config=task_config,
        )

        # pp.pprint(self.result.model_dump())
        # logger.debug(pp.pformat(self.result.model_dump()))
        self.computed = True

    def get_results(self, client: Optional["qcportal.FractalClient"] = None) -> AtomicResult:
        """Return results as Atomic-flavored QCSchema.

        NOTE: client removed (compared to psi4.driver.AtomicComputer)
        """
        if self.result:
            return self.result


class ManyBodyComputer(BaseComputerQCNG):

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
    embedding_charges: Optional[Dict[int, List[float]]] = Field(
        None,
        description="Atom-centered point charges to be used to speed up nbody-level convergence. Charges are placed on "
        "molecule fragments whose basis sets are not included in the computation. (An implication is that charges aren't "
        "invoked for bsse_type=cp.) Keys: 1-based index of fragment. Values: list of atom charges for that fragment.",
        # TODO: enforce point charge sum == fragment_charges val
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
        description=ManyBodyKeywords.__fields__["levels"].field_info.description
        + "Examples above are processed in the ManyBodyComputer, and once processed, only the values should be used. "
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
    task_list: Dict[str, Any] = {}  # MBETaskComputers] = {}
    qcmb_core: Optional[Any] = Field(
        None,
        description="Low-level interface",
    )

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
        return list(
            dict.fromkeys(
                [(BsseEnum[bt.lower()].value if bt.lower() in BsseEnum.__members__ else bt.lower()) for bt in v]
            )
        )

    # v2: @field_validator("embedding_charges")
    @validator("embedding_charges", pre=True)
    @classmethod
    # v2: def set_embedding_charges(cls, v: Any, info: FieldValidationInfo) -> Dict[int, List[float]]:
    def set_embedding_charges(cls, v, values):  # -> Dict[int, List[float]]:
        # print(f"hit embedding_charges validator with {v}", end="")
        nfr = len(values["molecule"].fragments)
        # v2: if len(v) != info.data["nfragments"]:
        if v and len(v) != nfr:
            raise ValueError(f"embedding_charges dict should have entries for each 1-indexed fragment ({nfr})")

        # print(f" ... setting embedding_charges={v}")
        return v

    # v2: @field_validator("return_total_data")
    @validator("return_total_data", always=True)
    @classmethod
    # v2: def set_return_total_data(cls, v: Optional[bool], info: FieldValidationInfo) -> bool:
    def set_return_total_data(cls, v: Optional[bool], values) -> bool:
        # print(f"hit return_total_data validator with {v}", end="")
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

        # print(f" ... setting rtd={rtd}")
        return rtd

    # v2: @field_validator("levels")
    @validator("levels", always=True)
    @classmethod
    # v2: def set_levels(cls, v: Any, info: FieldValidationInfo) -> Dict[Union[int, Literal["supersystem"]], str]:
    def set_levels(cls, v: Any, values) -> Dict[Union[int, Literal["supersystem"]], str]:
        # print(f"hit levels validator with {v}", end="")

        if v is None:
            pass
            # TODO levels = {plan.max_nbody: method}
            # v = {info.data["nfragments"]: "???method???"}
            # v = {len(info.data["molecule"].fragments): "???method???"}
            # v2: v = {len(info.data["molecule"].fragments): "(auto)"}
            v = {len(values["molecule"].fragments): "(auto)"}
        else:
            # rearrange bodies in order with supersystem last lest body count fail in organization loop below
            v = dict(sorted(v.items(), key=lambda item: 1000 if item[0] == "supersystem" else item[0]))

        # print(f" ... setting levels={v}")
        return v

    # TODO @computed_field(
    # TODO     description="Distribution of active n-body levels among model chemistry levels. All bodies in range "
    # TODO         "[1, self.max_nbody] must be present exactly once. Number of items in outer list is how many different "
    # TODO         "modelchems. Each inner list specifies what n-bodies to be run at the corresponding modelchem (e.g., "
    # TODO         "`[[1, 2]]` has max_nbody=2 and 1-body and 2-body contributions computed at the same level of theory; "
    # TODO         "`[[1], [2]]` has max_nbody=2 and 1-body and 2-body contributions computed at different levels of theory. "
    # TODO         "An entry 'supersystem' means all higher order n-body effects up to the number of fragments. The n-body "
    # TODO         "levels are effectively sorted in the outer list, and any 'supersystem' element is at the end.")
    # json_schema_extra={
    #    "shape": ["nmc", "<varies>"],
    # },
    @property
    def nbodies_per_mc_level(self) -> List[List[Union[int, Literal["supersystem"]]]]:
        # print(f"hit nbodies_per_mc_level", end="")

        # Organize nbody calculations into modelchem levels
        # * expand keys of `levels` into full lists of nbodies covered. save to plan, resetting max_nbody accordingly
        # * below, process values of `levels`, which are modelchem strings, into kwargs specs
        nbodies_per_mc_level = []
        prev_body = 0
        for nb in self.levels:
            nbodies = []
            if nb == "supersystem":
                nbodies.append(nb)
            elif nb != (prev_body + 1):
                for m in range(prev_body + 1, nb + 1):
                    nbodies.append(m)
            else:
                nbodies.append(nb)
            nbodies_per_mc_level.append(nbodies)
            prev_body = nb  # formerly buggy `+= 1`

        # print(f" ... setting nbodies_per_mc_level={nbodies_per_mc_level}")
        return nbodies_per_mc_level

    # v2: @field_validator("max_nbody")
    @validator("max_nbody", always=True)
    @classmethod
    # v2: def set_max_nbody(cls, v: Any, info: FieldValidationInfo) -> int:
    def set_max_nbody(cls, v: Any, values) -> int:
        # print(f"hit max_nbody validator with {v}", end="")
        # v2: levels_max_nbody = max(nb for nb in info.data["levels"] if nb != "supersystem")
        # v2: nfr = len(info.data["molecule"].fragments)
        levels_max_nbody = max(nb for nb in values["levels"] if nb != "supersystem")
        nfr = len(values["molecule"].fragments)
        # print(f" {levels_max_nbody=} {nfr=}", end="")

        if len(set(values["levels"].values())) != len(values["levels"]):
            raise ValueError("Cannot have duplicate model chemistries in levels.")

        # ALT if v == -1:
        if v is None:
            v = levels_max_nbody
        elif v < 0 or v > nfr:
            raise ValueError(f"max_nbody={v} should be between 1 and {nfr}.")
        elif v != levels_max_nbody:
            # raise ValueError(f"levels={levels_max_nbody} contradicts user max_nbody={v}.")
            # TODO reconsider logic. move this from levels to here?
            # v2: info.data["levels"] = {v: "(auto)"}
            values["levels"] = {v: "(auto)"}
        else:
            pass
            # TODO once was           return min(v, nfragments)

        # print(f" ... setting max_nbody={v}")
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
        # print(f"hit supersystem_ie_only validator with {v}", end="")
        sio = v
        # v2: _nfr = len(info.data["molecule"].fragments)
        _nfr = len(values["molecule"].fragments)

        # v2: _max_nbody = info.data["max_nbody"]
        # get(..., None) b/c in v1, all fields processed even if max_nbody previously failed
        _max_nbody = values.get("max_nbody", None)
        if (sio is True) and (_max_nbody != _nfr):
            raise ValueError(f"Cannot skip intermediate n-body jobs when max_nbody={_max_nbody} != nfragments={_nfr}.")

        if (sio is True) and ("vmfc" in values["bsse_type"]):
            raise ValueError(
                f"Cannot skip intermediate n-body jobs when VMFC in bsse_type={values['bsse_type']}. Use CP instead."
            )

        # print(f" ... setting {sio=}")
        return sio

    @classmethod
    def from_manybodyinput(cls, input_model: ManyBodyInput, build_tasks: bool = True):

        computer_model = cls(
            molecule=input_model.molecule,
            driver=input_model.specification.driver,
            # v2: **input_model.specification.keywords.model_dump(exclude_unset=True),
            **input_model.specification.keywords.dict(exclude_unset=True),
            input_data=input_model,  # storage, to reconstitute ManyBodyResult
        )
        nb_per_mc = computer_model.nbodies_per_mc_level

        # print("\n<<<  (ZZ 1) QCEngine harness ManyBodyComputerQCNG.from_qcschema_ben  >>>")
        # v2: pprint.pprint(computer_model.model_dump(), width=200)
        # pprint.pprint(computer_model.dict(), width=200)
        # print(f"nbodies_per_mc_level={nb_per_mc}")

        comp_levels = {}
        for mc_level_idx, mtd in enumerate(computer_model.levels.values()):
            for lvl1 in nb_per_mc[mc_level_idx]:
                key = "supersystem" if lvl1 == "supersystem" else int(lvl1)
                comp_levels[key] = mtd

        specifications = {}
        for mtd, spec in computer_model.input_data.specification.specification.items():
            spec = spec.dict()
            specifications[mtd] = {}
            specifications[mtd]["program"] = spec.pop("program")
            specifications[mtd]["specification"] = spec
            specifications[mtd]["specification"][
                "driver"
            ] = computer_model.driver  # overrides atomic driver with mb driver
            specifications[mtd]["specification"].pop("schema_name", None)

        computer_model.qcmb_core = ManyBodyCore(
            computer_model.molecule,
            computer_model.bsse_type,
            comp_levels,
            return_total_data=computer_model.return_total_data,
            supersystem_ie_only=computer_model.supersystem_ie_only,
            embedding_charges=computer_model.embedding_charges,
        )

        # check that core and computer storage are consistent in mc ordering and grouping and nbody levels
        assert (
            list(computer_model.qcmb_core.nbodies_per_mc_level.values()) == computer_model.nbodies_per_mc_level
        ), f"CORE {computer_model.qcmb_core.nbodies_per_mc_level.values()} != COMPUTER {computer_model.nbodies_per_mc_level}"
        assert list(computer_model.qcmb_core.nbodies_per_mc_level.keys()) == list(
            computer_model.levels.values()
        ), f"CORE {computer_model.qcmb_core.nbodies_per_mc_level.keys()} != COMPUTER {computer_model.levels.values()}"

        if not build_tasks:
            return computer_model

        try:
            import qcengine as qcng
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "Python module qcengine not found. Solve by installing it: "
                "`conda install qcengine -c conda-forge` or `pip install qcengine`"
            )

        component_properties = {}
        component_results = {}

        for chem, label, imol in computer_model.qcmb_core.iterate_molecules():
            inp = AtomicInput(molecule=imol, **specifications[chem]["specification"])
            # inp = AtomicInput(molecule=imol, **specifications[chem]["specification"], extras={"psiapi": True})  # faster for p4

            if imol.extras.get("embedding_charges"):  # or test on self.embedding_charges ?
                if specifications[chem]["program"] == "psi4":
                    charges = imol.extras["embedding_charges"]
                    fkw = inp.keywords.get("function_kwargs", {})
                    fkw.update({"external_potentials": charges})
                    inp.keywords["function_kwargs"] = fkw
                else:
                    raise RuntimeError(
                        f"Don't know how to handle external charges in {specifications[chem]['program']}"
                    )

            _, real, bas = delabeler(label)
            result = qcng.compute(inp, specifications[chem]["program"])
            component_results[label] = result

            if not result.success:
                # print(result.error.error_message)
                raise RuntimeError("Calculation did not succeed! Error:\n" + result.error.error_message)

            # pull out stuff
            props = {"energy", "gradient", "hessian"}

            component_properties[label] = {}

            for p in props:
                if hasattr(result.properties, f"return_{p}"):
                    v = getattr(result.properties, f"return_{p}")
                    # print(f"  {label} {p}: {v}")
                    if v is not None:
                        component_properties[label][p] = v

        # print("\n<<<  (ZZ 2) QCEngine harness ManyBodyComputerQCNG.from_qcschema_ben component_properties  >>>")
        # with np.printoptions(precision=6, suppress=True):
        #     pprint.pprint(component_properties, width=200)

        analyze_back = computer_model.qcmb_core.analyze(component_properties)
        analyze_back["nbody_number"] = len(component_properties)
        # print("\n<<<  (ZZ 3) QCEngine harness ManyBodyComputerQCNG.from_qcschema_ben analyze_back  >>>")
        # pprint.pprint(analyze_back, width=200)

        return computer_model.get_results(external_results=analyze_back, component_results=component_results)

    def plan(self):
        # uncalled function
        return [t.plan() for t in self.task_list.values()]

    def compute(self, client: Optional["qcportal.FractalClient"] = None) -> None:
        """Run quantum chemistry.

        NOTE: client logic removed (compared to psi4.driver.ManyBodyComputer)
        """
        for t in self.task_list.values():
            t.compute(client=client)

    def get_results(
        self, external_results: Dict, component_results: Dict, client: Optional["qcportal.FractalClient"] = None
    ) -> ManyBodyResult:
        """Return results as ManyBody-flavored QCSchema."""

        ret_energy = external_results.pop("ret_energy")
        ret_ptype = ret_energy if self.driver == "energy" else external_results.pop(f"ret_{self.driver.name}")
        ret_gradient = external_results.pop("ret_gradient", None)
        nbody_number = external_results.pop("nbody_number")
        component_properties = external_results.pop("component_properties")
        stdout = external_results.pop("stdout")

        properties = {
            "calcinfo_nmc": len(self.nbodies_per_mc_level),
            "calcinfo_nfr": self.nfragments,  # or len(self.molecule.fragments)
            "calcinfo_natom": len(self.molecule.symbols),
            "calcinfo_nmbe": nbody_number,
            "nuclear_repulsion_energy": self.molecule.nuclear_repulsion_energy(),
            "return_energy": ret_energy,
        }

        if self.driver == "gradient":
            properties["return_gradient"] = ret_ptype
        elif self.driver == "hessian":
            properties["return_gradient"] = ret_gradient
            properties["return_hessian"] = ret_ptype

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

        # print("QCVARS PRESCREEN")
        # pp.pprint(qcvars)

        # v2: component_results = self.model_dump()['task_list']  # TODO when/where include the indiv outputs
        # ?component_results = self.dict()['task_list']  # TODO when/where include the indiv outputs
        #        for k, val in component_results.items():
        #            val['molecule'] = val['molecule'].to_schema(dtype=2)

        nbody_model = ManyBodyResult(
            **{
                "input_data": self.input_data,
                #'molecule': self.molecule,
                # v2: 'properties': {**atprop.model_dump(), **properties},
                "properties": {**external_results["results"], **properties},
                "component_properties": component_properties,
                "component_results": component_results,
                "provenance": provenance_stamp(__name__),
                "return_result": ret_ptype,
                "stdout": stdout,
                "success": True,
            }
        )

        #        logger.debug('\nNBODY QCSchema:\n' + pp.pformat(nbody_model.model_dump()))

        return nbody_model
