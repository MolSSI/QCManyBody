from __future__ import annotations

from enum import Enum, IntEnum
from typing import Any, Dict, List, Optional, Literal, Tuple, Union

# v2: from pydantic import create_model, Field, field_validator, FieldValidationInfo
try:
    from pydantic.v1 import create_model, Field, validator
except ImportError:
    from pydantic import create_model, Field, validator

from qcelemental.models.types import Array
#from .basemodels import ExtendedConfigDict, ProtoModel
from qcelemental.models.common_models import Model
from qcelemental.models.molecule import Molecule
from qcelemental.models.results import AtomicResultProtocols
from qcelemental.models import DriverEnum, ProtoModel, Provenance


# ====  Misplaced & Next Models  ================================================

class AtomicSpecification(ProtoModel):
    """Specification for a single point QC calculation"""

    keywords: Dict[str, Any] = Field({}, description="The program specific keywords to be used.")
    program: str = Field(..., description="The program for which the Specification is intended.")

    schema_name: Literal["qcschema_atomicspecification"] = "qcschema_atomicspecification"
    driver: DriverEnum = Field(..., description=DriverEnum.__doc__)
    model: Model = Field(..., description=Model.__doc__)
    protocols: AtomicResultProtocols = Field(
        AtomicResultProtocols(),
        description=AtomicResultProtocols.__doc__,
    )


class ResultBase(ProtoModel):
    """Base class for all result classes"""

    #input_data: InputBase = Field(..., description=InputBase.__doc__)
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


# ====  Inputs  =================================================================

class BsseEnum(str, Enum):
#class BsseEnum(IntEnum):
    """Available basis-set superposition error (BSSE) treatments."""

    nocp = "nocp"  # plain supramolecular interaction energy
    cp = "cp"      # Boys-Bernardi counterpoise correction; site-site functional counterpoise (SSFC)
    vmfc = "vmfc"  # Valiron-Mayer function counterpoise
    ssfc = "cp"


FragBasIndex = Tuple[Tuple[int], Tuple[int]]


class ManyBodyKeywords(ProtoModel):
    """The many-body-specific keywords for user control."""

    bsse_type: List[BsseEnum] = Field(
        [BsseEnum.cp],
        # definitive description
        description="Requested BSSE treatments. First in list determines which interaction or total "
            "energy/gradient/Hessian returned.",
    )
    embedding_charges: Dict[int, List[float]] = Field(
        {},
        description="Atom-centered point charges to be used on molecule fragments whose basis sets are not included in "
            "the computation. Keys: 1-based index of fragment. Values: list of atom charges for that fragment.",
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
            "effects up to the number of fragments. A method fills in for any lower unlisted nbody levels. Note that if "
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
            "in the MBE up through ``max_nbody``. When True (only allowed for ``max_nbody = nfragments``), only compute "
            "enough for the overall interaction/total energy: max_nbody-body and 1-body. When True, properties "
            "``INTERACTION {driver} THROUGH {max_nbody}-BODY`` will always be available; ``TOTAL {driver} THROUGH "
            "{max_nbody}-BODY`` will be available depending on ``return_total_data``; and ``{max_nbody}-BODY "
            "CONTRIBUTION TO {driver}`` won't be available (except for dimers). This keyword produces no savings for a "
            "two-fragment molecule. But for the interaction energy of a three-fragment molecule, for example, 2-body "
            "subsystems can be skipped with ``supersystem_ie_only=True`` Do not use with ``vmfc`` in ``bsse_type``"
            "as it cannot produce savings."
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
        return list(dict.fromkeys([(BsseEnum[bt.lower()].value if bt.lower() in BsseEnum.__members__ else bt.lower()) for bt in v]))


class ManyBodySpecification(ProtoModel):

    schema_name: Literal["qcschema_manybodyspecification"] = "qcschema_manybodyspecification"
    #provenance: Provenance = Field(Provenance(**provenance_stamp(__name__)), description=Provenance.__doc__)
    keywords: ManyBodyKeywords = Field(..., description=ManyBodyKeywords.__doc__)
    #program: str = Field(..., description="The program for which the Specification is intended.")
    driver: DriverEnum = Field(
        ...,
        description="The computation driver; i.e., energy, gradient, hessian.",
    )
    #specification: Union[AtomicSpecification, Dict[str, AtomicSpecification]] = Field(
    specification: Dict[str, AtomicSpecification] = Field(
        ...,
        description="??? TODO expand to cbs, fd",
    )

    # v2: @field_validator("specification", mode="before")
    @validator("specification", pre=True)
    @classmethod
    def set_specification(cls, v: Any) -> Dict[str, AtomicSpecification]:
        #print(f"hit atomicspecification validator with {type(v)=} {v}", end="")
        # v could be model instance or dict
        if isinstance(v, AtomicSpecification) or "model" in v:
            v = {"(auto)": v}
        #print(f" ... setting v={v}")
        return v


class ManyBodyInput(ProtoModel):

    schema_name: Literal["qcschema_manybodyinput"] = "qcschema_manybodyinput"
    #provenance: Provenance = Field(Provenance(**provenance_stamp(__name__)), description=Provenance.__doc__)
    specification: ManyBodySpecification = Field(
        ...,
        description="???",
    )
    molecule: Molecule = Field(
        ...,
        description="Target molecule for many-body expansion (MBE) or interaction energy (IE) analysis.",
    )
    #protocols


# ====  Protocols  ==============================================================

class ManyBodyProtocolEnum(str, Enum):
    """
    Which atomic evaluations to keep in a many body evaluation.
    """

    all = "all"
    all_real = "all_real"
    largest_body = "largest_body"
    none = "none"


class ManyBodyProtocols(ProtoModel):
    """
    Protocols regarding the manipulation of a ManyBody output data.
    """

    atomics: ManyBodyProtocolEnum = Field(
        ManyBodyProtocolEnum.all,
        description=str(ManyBodyProtocolEnum.__doc__),
    )

    # v2: model_config = ExtendedConfigDict(force_skip_defaults=True)
    class Config:
        force_skip_defaults = True


# ====  Properties  =============================================================

# class ManyBodyResultProperties defined through create_model

manybodyresultproperties_doc = """
    Named properties of manybody computations following the MolSSI QCSchema.

    All arrays are stored flat but must be reshapable into the dimensions in attribute ``shape``, with abbreviations as follows:

    * nat: number of atoms = :attr:`~qcelemental.models.ManyBodyResultProperties.calcinfo_natom`
    * nmc: number of model chemistries = :attr:`~qcelemental.models.ManyBodyResultProperties.calcinfo_nmc`
    * nfr: number of fragments = :attr:`~qcelemental.models.ManyBodyResultProperties.calcinfo_nfr`
    * nmbe: number of jobs = :attr:`~qcelemental.models.ManyBodyResultProperties.calcinfo_nmbe`
    """

MAX_NBODY = 5  # 5 covers tetramers

json_schema_extras = {
    "energy": {"units": "E_h"},
    "gradient": {"units": "E_h/a0", "shape": ["nat", 3]},
    "Hessian": {"units": "E_h/a0^2", "shape": ["nat" * 3, "nat" * 3]},
}

mbprop = {}

# ========  Calcinfo  ===========================================================

mbprop["calcinfo_nmc"] = (
        Optional[int],
        Field(
            None,
            description="The number of model chemistries applied to n-body levels of the computation.",
        ))

mbprop["calcinfo_nfr"] = (
        Optional[int],
        Field(
            None,
            description="The number of fragments in the molecule for the computation.",
        ))

mbprop["calcinfo_natom"] = (
        Optional[int],
        Field(
            None,
            description="The number of atoms in the computation.",
        ))  # alias nat

mbprop["calcinfo_nmbe"] = (
        Optional[int],
        Field(
            None,
            description="The number of real/ghost molecule patterns for the computation.",
        ))  # alias NBODY NUMBER

# ========  Canonical  ==========================================================

mbprop["nuclear_repulsion_energy"] = (
        Optional[float],
        Field(
            None,
            description="The nuclear repulsion energy.",
        ))

# ret_energy
mbprop["return_energy"] = (
        Optional[float],
        Field(
            None,
            description=f"The interaction energy of the requested method: IE or total (depending on return_total_data) with cp/nocp/vmfc treatment (dep. on first of bsse_type). Always available. Identical to :attr:`~qcelemental.models.ManyBodyResult.return_result` for :attr:`~qcelemental.models.AtomicInput.driver`\\ =\\ :attr:`~qcelemental.models.DriverEnum.energy` computations.",
            json_schema_extra={"units": "E_h"},
        ))

# ret_gradient
mbprop["return_gradient"] = (
        Optional[Array[float]],
        Field(
            None,
            description=f"The interaction gradient of the requested method: IE or total (depending on return_total_data) with cp/nocp/vmfc treatment (dep. on first of bsse_type). Available when driver is g/h. Identical to :attr:`~qcelemental.models.ManyBodyResult.return_result` for :attr:`~qcelemental.models.AtomicInput.driver`\\ =\\ :attr:`~qcelemental.models.DriverEnum.gradient` computations.",
            json_schema_extra=json_schema_extras["gradient"],
        ))

# ret_hessian
mbprop["return_hessian"] = (
        Optional[Array[float]],
        Field(
            None,
            description=f"The interaction Hessian of the requested method: IE or total (depending on return_total_data) with cp/nocp/vmfc treatment (dep. on first of bsse_type). Available when driver is h. Identical to :attr:`~qcelemental.models.ManyBodyResult.return_result` for :attr:`~qcelemental.models.AtomicInput.driver`\\ =\\ :attr:`~qcelemental.models.DriverEnum.hessian` computations.",
            json_schema_extra=json_schema_extras["Hessian"],
        ))

# ========  CP E/G/H summary data  ==============================================

for singular in ["energy", "gradient", "Hessian"]:
    jse = json_schema_extras[singular]
    egh = singular.lower()
    typ = float if singular == "energy" else Array[float]
    plural = "energies" if singular == "energy" else singular + "s"
    availability_of_derivative = {
        "energy": "",
        "gradient": " & driver is g/h",
        "Hessian": " & driver is h",
    }[singular]

    # CP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"cp_corrected_total_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"MBE sum of subsystems of {nb}-body or fewer (cumulative); summed are total {plural} w/ cp treatment. Available when cp in bsse_type & rtd=T & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # CP-CORRECTED TOTAL ENERGY
    mbprop[f"cp_corrected_total_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"Best available total {singular} with cp treatment: cp_corrected_total_{egh}_through_{{max_nbody}}_body. Available when cp in bsse_type & rtd=T{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # CP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"cp_corrected_interaction_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less 1-body total data for cumulative IE; inputs are total {plural} with cp treatment. Available when when cp in bsse_type & max_nbody>={nb}{availability_of_derivative}. The 1-body quantity is zero by definition.",
                json_schema_extra=jse,
            ))

    # CP-CORRECTED INTERACTION ENERGY
    mbprop[f"cp_corrected_interaction_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"Best available interaction {singular} with cp treatment: cp_corrected_interaction_{egh}_through_{{max_nbody}}_body. Available when cp in bsse_type{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # CP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY
    for nb in range(2, MAX_NBODY):
        mbprop[f"cp_corrected_{nb}_body_contribution_to_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less ({nb}-1)-body data for partial IE; inputs are total {plural} w/ cp treat. Available when cp in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

# ========  NOCP E/G/H summary data  ============================================

    # NOCP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"nocp_corrected_total_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"MBE sum of subsystems of {nb}-body or fewer (cumulative); summed are total {plural} without cp treatment. Available when nocp in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # NOCP-CORRECTED TOTAL ENERGY
    mbprop[f"nocp_corrected_total_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"Best available total {singular} without cp treatment: nocp_corrected_total_{egh}_through_{{max_nbody}}_body. Available when nocp in bsse_type{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # NOCP-CORRECTED INTERATION ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"nocp_corrected_interaction_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less 1-body total data for cumulative IE; inputs are total {plural} without cp treatment. Available when when nocp in bsse_type & max_nbody>={nb}{availability_of_derivative}. The 1-body quantity is zero by definition.",
                json_schema_extra=jse,
            ))

# TODO note htat TOT 1BODY cp=nocp=vmfc
# TODO note that summ INTERACTION ENERGY props (w/o explicit -BODY) return 0.0 for max_nbody=1 for completeness

    # NOCP-CORRECTED INTERACTION ENERGY
    mbprop[f"nocp_corrected_interaction_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"Best available interaction {singular} without cp treatment: nocp_corrected_interaction_{egh}_through_{{max_nbody}}_body. Available when nocp in bsse_type{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # NOCP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY
    for nb in range(2, MAX_NBODY):
        mbprop[f"nocp_corrected_{nb}_body_contribution_to_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less ({nb}-1)-body data for partial IE; inputs are total {plural} w/o cp treatment. Available when nocp in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

# ========  VMFC E/G/H summary data  ============================================

    # VMFC-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"vmfc_corrected_total_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"MBE sum of subsystems of {nb}-body or fewer (cumulative); summed are total {plural} with vmfc treatment. Available when vmfc in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # VMFC-CORRECTED TOTAL ENERGY
    mbprop[f"vmfc_corrected_total_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"Best available total {singular} with vmfc treatment: vmfc_corrected_total_{egh}_through_{{max_nbody}}_body. Available when vmfc in bsse_type{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # VMFC-CORRECTED INTERATION ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"vmfc_corrected_interaction_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less 1-body total data for cumulative IE; inputs are total {plural} w/ vmfc treatment. Available when when vmfc in bsse_type & max_nbody>={nb}{availability_of_derivative}. The 1-body quantity is zero by definition.",
                json_schema_extra=jse,
            ))

    # VMFC-CORRECTED INTERACTION ENERGY
    mbprop[f"vmfc_corrected_interaction_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"Best available interaction {singular} with vmfc treatment: vmfc_corrected_interaction_{egh}_through_{{max_nbody}}_body. Available when vmfc in bsse_type{availability_of_derivative}.",
                json_schema_extra=jse,
            ))

    # VMFC-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY
    for nb in range(2, MAX_NBODY):
        mbprop[f"vmfc_corrected_{nb}_body_contribution_to_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less ({nb}-1)-body total data for partial IE; inputs are total {plural} w/ vmfc treatment. Available when vmfc in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ))


ManyBodyResultProperties = create_model(
    "ManyBodyResultProperties",
    #__config__=ConfigDict(title='abc'),
    #__doc__=manybodyresultproperties_doc,  # needs later pydantic
    __base__=ProtoModel,
    **mbprop,
)


# ====  Results  ================================================================

# ManyBodyResult(ManyBodyInput):
class ManyBodyResult(SuccessfulResultBase):

    schema_name: Literal["qcschema_manybodyresult"] = "qcschema_manybodyresult"
    schema_version: Literal[1] = Field(
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    id: Optional[str] = Field(None, description="The optional ID for the object.")
    extras: Dict[str, Any] = Field(
        {},
        description="Additional information to bundle with the object. Use for schema development and scratch space.",
    )

    provenance: Provenance = Field(..., description=str(Provenance.__doc__))
    input_data: ManyBodyInput = Field(
        ...,
    )
    success: bool = Field(
        ...,
        description="A boolean indicator that the operation succeeded or failed. Allows programmatic assessment of "
        "all results regardless of if they failed or succeeded by checking `result.success`.",
    )
    properties: ManyBodyResultProperties = Field(..., description=str(ManyBodyResultProperties.__doc__))
    return_result: Union[float, Array[float], Dict[str, Any]] = Field(
        ...,
        description="The primary return specified by the :attr:`~qcelemental.models.AtomicInput.driver` field. Scalar if energy; array if gradient or hessian; dictionary with property keys if properties.",
    )
    stdout: Optional[str] = Field(
        None,
        description="The primary logging output of the program, whether natively standard output or a file. Presence vs. absence (or null-ness?) configurable by protocol.",
    )
    stderr: Optional[str] = Field(None, description="The standard error of the program execution.")
    success: Literal[True] = Field(True, description="Always `True` for a successful result")
