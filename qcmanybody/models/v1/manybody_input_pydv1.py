from __future__ import annotations

from enum import Enum, IntEnum
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional, Tuple, Union

# v2: from pydantic import create_model, Field, field_validator, FieldValidationInfo
try:
    from pydantic.v1 import Field, create_model, validator
except ImportError:
    from pydantic import create_model, Field, validator

from qcelemental.models import DriverEnum, ProtoModel, Provenance

# from .basemodels import ExtendedConfigDict, ProtoModel
from qcelemental.models.common_models import Model
from qcelemental.models.molecule import Molecule
from qcelemental.models.results import AtomicResultProperties, AtomicResultProtocols
from qcelemental.models.types import Array

# ====  Misplaced & Next Models  ================================================


class AtomicSpecification(ProtoModel):
    """Specification for a single point QC calculation"""

    keywords: Dict[str, Any] = Field({}, description="The program specific keywords to be used.")
    program: str = Field(..., description="The program for which the Specification is intended.")

    schema_name: Literal["qcschema_atomicspecification"] = "qcschema_atomicspecification"
    schema_version: Literal[1] = Field(
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )

    driver: DriverEnum = Field(..., description=DriverEnum.__doc__)
    model: Model = Field(..., description=Model.__doc__)
    protocols: AtomicResultProtocols = Field(
        AtomicResultProtocols(),
        description=AtomicResultProtocols.__doc__,
    )
    extras: Dict[str, Any] = Field(
        {},
        description="Additional information to bundle with the computation. Use for schema development and scratch space.",
    )


class ResultBase(ProtoModel):
    """Base class for all result classes"""

    # input_data: InputBase = Field(..., description=InputBase.__doc__)
    input_data: Any
    success: bool = Field(
        ...,
        description="A boolean indicator that the operation succeeded or failed. Allows programmatic assessment of "
        "all results regardless of if they failed or succeeded by checking `result.success`.",
    )
    stdout: Optional[str] = Field(
        None,
        description="The primary logging output of the program, whether natively standard output or a file. Presence vs. absence (or null-ness?) configurable by protocol.",
    )
    stderr: Optional[str] = Field(None, description="The standard error of the program execution.")


class SuccessfulResultBase(ResultBase):
    """Base object for any successful result"""

    success: Literal[True] = Field(True, description="Always `True` for a successful result")


# ====  Protocols  ==============================================================


class ComponentResultsProtocolEnum(str, Enum):
    r"""Which component results to preserve in a many body result; usually AtomicResults."""

    all = "all"
    # max_nbody = "max_nbody"
    none = "none"


class ManyBodyProtocols(ProtoModel):
    """
    Protocols regarding the manipulation of a ManyBody output data.
    """

    component_results: ComponentResultsProtocolEnum = Field(
        ComponentResultsProtocolEnum.none, description=str(ComponentResultsProtocolEnum.__doc__)
    )

    # v2: model_config = ExtendedConfigDict(force_skip_defaults=True)
    class Config:
        force_skip_defaults = True


# ====  Inputs  =================================================================


class BsseEnum(str, Enum):
    """Available basis-set superposition error (BSSE) treatments."""

    nocp = "nocp"  # plain supramolecular interaction energy
    cp = "cp"  # Boys-Bernardi counterpoise correction; site-site functional counterpoise (SSFC)
    vmfc = "vmfc"  # Valiron-Mayer function counterpoise
    ssfc = "cp"
    mbe = "nocp"
    none = "nocp"

    def formal(self):
        return {
            "nocp": "Non-Counterpoise Corrected",
            "cp": "Counterpoise Corrected",
            "vmfc": "Valiron-Mayer Function Counterpoise",
        }[self]

    def abbr(self):
        return {
            "nocp": "NoCP",
            "cp": "CP",
            "vmfc": "VMFC",
        }[self]


FragBasIndex = Tuple[Tuple[int], Tuple[int]]


class ManyBodyKeywords(ProtoModel):
    """The many-body-specific keywords for user control."""

    schema_name: Literal["qcschema_manybodykeywords"] = "qcschema_manybodykeywords"
    schema_version: Literal[1] = Field(
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    bsse_type: List[BsseEnum] = Field(
        [BsseEnum.cp],
        # definitive description
        description="Requested BSSE treatments. First in list determines which interaction or total "
        "energy/gradient/Hessian returned.",
    )
    embedding_charges: Optional[Dict[int, List[float]]] = Field(
        None,
        description="Atom-centered point charges to be used on molecule fragments whose basis sets are not included in "
        "the computation. Keys: 1-based index of fragment. Values: list of atom charges for that fragment. "
        "At present, QCManyBody will only accept non-None values of this keyword if environment variable "
        "QCMANYBODY_EMBEDDING_CHARGES is set.",
        # TODO embedding charges should sum to fragment charge, right? enforce?
        # TODO embedding charges irrelevant to CP (basis sets always present)?
        json_schema_extra={
            "shape": ["nfr", "<varies: nat in ifr>"],
        },
    )
    return_total_data: Optional[bool] = Field(
        None,
        validate_default=True,
        # definitive description
        description="When True, returns the total data (energy/gradient/Hessian) of the system, otherwise returns "
        "interaction data. Default is False for energies, True for gradients and Hessians. Note that the calculation "
        "of counterpoise corrected total energies implies the calculation of the energies of monomers in the monomer "
        "basis, hence specifying ``return_total_data = True`` may carry out more computations than "
        "``return_total_data = False``. For gradients and Hessians, ``return_total_data = False`` is rarely useful.",
    )
    levels: Optional[Dict[Union[int, Literal["supersystem"]], str]] = Field(
        None,
        # definitive description. appended in Computer
        description="Dictionary of different levels of theory for different levels of expansion. Note that the primary "
        "method_string is not used when this keyword is given. ``supersystem`` computes all higher order n-body "
        "effects up to the number of fragments; this higher-order correction uses the nocp basis, regardless of "
        "bsse_type. A method fills in for any lower unlisted nbody levels. Note that if "
        "both this and max_nbody are provided, they must be consistent. Examples: "
        "SUPERSYSTEM definition suspect"
        "* {1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'} "
        "* {2: 'ccsd(t)/cc-pvdz', 3: 'mp2'} "
        "* Now invalid: {1: 2, 2: 'ccsd(t)/cc-pvdz', 3: 'mp2'} ",
    )
    max_nbody: Optional[int] = Field(
        None,
        validate_default=True,
        # definitive description
        description="Maximum number of bodies to include in the many-body treatment. Possible: max_nbody <= nfragments. "
        "Default: max_nbody = nfragments.",
    )
    supersystem_ie_only: Optional[bool] = Field(
        False,
        validate_default=True,
        # definitive description
        description="Target the supersystem total/interaction energy (IE) data over the many-body expansion (MBE) "
        "analysis, thereby omitting intermediate-body calculations. When False (default), compute each n-body level "
        "in the MBE up through ``max_nbody``. When True (only allowed for ``max_nbody = nfragments`` ), only compute "
        "enough for the overall interaction/total energy: max_nbody-body and 1-body. When True, properties "
        "``INTERACTION {driver} THROUGH {max_nbody}-BODY`` will always be available; "
        "``TOTAL {driver} THROUGH {max_nbody}-BODY`` will be available depending on ``return_total_data`` ; and "
        "``{max_nbody}-BODY CONTRIBUTION TO {driver}`` won't be available (except for dimers). This keyword produces "
        "no savings for a two-fragment molecule. But for the interaction energy of a three-fragment molecule, for "
        "example, 2-body subsystems can be skipped with ``supersystem_ie_only=True``. Do not use with ``vmfc`` in "
        "``bsse_type`` as it cannot produce savings.",
    )

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


class ManyBodySpecification(ProtoModel):
    """Combining the what (ManyBodyKeywords) with the how (AtomicSpecification)."""

    schema_name: Literal["qcschema_manybodyspecification"] = "qcschema_manybodyspecification"
    schema_version: Literal[1] = Field(
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    # provenance: Provenance = Field(Provenance(**provenance_stamp(__name__)), description=Provenance.__doc__)
    keywords: ManyBodyKeywords = Field(..., description=ManyBodyKeywords.__doc__)
    # program: str = Field(..., description="The program for which the Specification is intended.")  # TODO is qcmanybody
    protocols: ManyBodyProtocols = Field(ManyBodyProtocols(), description=str(ManyBodyProtocols.__doc__))
    driver: DriverEnum = Field(
        ...,
        description="The computation driver; i.e., energy, gradient, hessian.",
    )
    # specification: Union[AtomicSpecification, Dict[str, AtomicSpecification]] = Field(
    specification: Dict[str, AtomicSpecification] = Field(
        ...,
        description="??? TODO expand to cbs, fd",
    )
    extras: Dict[str, Any] = Field(
        {},
        description="Additional information to bundle with the computation. Use for schema development and scratch space.",
    )

    # v2: @field_validator("specification", mode="before")
    @validator("specification", pre=True)
    @classmethod
    def set_specification(cls, v: Any) -> Dict[str, AtomicSpecification]:
        # print(f"hit atomicspecification validator with {type(v)=} {v}", end="")
        # v could be model instance or dict
        if isinstance(v, AtomicSpecification) or "model" in v:
            v = {"(auto)": v}
        # print(f" ... setting v={v}")
        return v


class ManyBodyInput(ProtoModel):
    """Combining the what and how (ManyBodySpecification) with the who (Molecule)."""

    schema_name: Literal["qcschema_manybodyinput"] = "qcschema_manybodyinput"
    schema_version: Literal[1] = Field(
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    # provenance: Provenance = Field(Provenance(**provenance_stamp(__name__)), description=Provenance.__doc__)
    specification: ManyBodySpecification = Field(
        ...,
        description="???",
    )
    molecule: Molecule = Field(
        ...,
        description="Target molecule for many-body expansion (MBE) or interaction energy (IE) analysis.",
    )
    extras: Dict[str, Any] = Field(
        {},
        description="Additional information to bundle with the computation. Use for schema development and scratch space.",
    )
