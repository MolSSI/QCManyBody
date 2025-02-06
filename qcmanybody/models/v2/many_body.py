from __future__ import annotations

import os
import re
from enum import Enum, IntEnum
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional, Tuple, Union

try:
    from pydantic import ValidationInfo
except ImportError:
    from pydantic import FieldValidationInfo as ValidationInfo

from pydantic import Field, create_model, field_validator, model_validator
from qcelemental.models.v2 import (  # Array,
    AtomicProperties,
    AtomicProtocols,
    AtomicResult,
    AtomicSpecification,
    DriverEnum,
    Model,
    Molecule,
    ProtoModel,
    Provenance,
)
from qcelemental.models.v2.basemodels import ExtendedConfigDict, ProtoModel, check_convertible_version
from qcelemental.models.v2.types import Array  # return to above once qcel corrected

from ...utils import provenance_stamp

# ====  Misplaced & Next Models  ================================================


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


class ClusterResultsProtocolEnum(str, Enum):
    r"""Which component results to preserve in a many body result; usually AtomicResults."""

    all = "all"
    # max_nbody = "max_nbody"
    none = "none"


class ManyBodyProtocols(ProtoModel):
    """
    Protocols regarding the manipulation of a ManyBody output data.
    """

    schema_name: Literal["qcschema_many_body_protocols"] = "qcschema_many_body_protocols"

    cluster_results: ClusterResultsProtocolEnum = Field(
        ClusterResultsProtocolEnum.none, description=str(ClusterResultsProtocolEnum.__doc__)
    )

    model_config = ExtendedConfigDict(force_skip_defaults=True)

    def convert_v(
        self, target_version: int, /
    ) -> Union["qcmanybody.models.v1.ManyBodyProtocols", "qcmanybody.models.v2.ManyBodyProtocols"]:
        """Convert to instance of particular QCSchema version."""
        import qcmanybody as qcmb

        if check_convertible_version(target_version, error="ManyBodyProtocols") == "self":
            return self

        dself = self.model_dump()
        if target_version == 1:
            # serialization is compact, so use model to assure value
            dself.pop("cluster_results", None)
            dself["component_results"] = self.cluster_results.value

            self_vN = qcmb.models.v1.ManyBodyProtocols(**dself)
        else:
            assert False, target_version

        return self_vN


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

    schema_name: Literal["qcschema_many_body_keywords"] = "qcschema_many_body_keywords"

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

    @field_validator("bsse_type", mode="before")
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

    schema_name: Literal["qcschema_many_body_specification"] = "qcschema_many_body_specification"

    keywords: ManyBodyKeywords = Field(..., description=ManyBodyKeywords.__doc__)
    program: str = Field(
        "", description="Many Body Expansion CMS code / QCEngine procedure with which to run the MB decomposition."
    )
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

    @field_validator("specification", mode="before")
    @classmethod
    def set_specification(cls, v: Any) -> Dict[str, AtomicSpecification]:
        # print(f"hit atomicspecification validator with {type(v)=} {v}", end="")
        # v could be model instance or dict
        if isinstance(v, AtomicSpecification) or "model" in v:
            v = {"(auto)": v}
        # print(f" ... setting v={v}")
        return v

    @field_validator("program")
    @classmethod
    def _check_procedure(cls, v):
        return v.lower()

    def convert_v(
        self, target_version: int, /
    ) -> Union["qcmanybody.models.v1.ManyBodySpecification", "qcmanybody.models.v2.ManyBodySpecification"]:
        """Convert to instance of particular QCSchema version."""
        import qcmanybody as qcmb

        if check_convertible_version(target_version, error="ManyBodySpecification") == "self":
            return self

        dself = self.model_dump()
        if target_version == 1:
            dself.pop("schema_name")
            dself["keywords"].pop("schema_name")
            try:
                dself["specification"].pop("schema_name")
            except KeyError:
                for spec in dself["specification"].values():
                    spec.pop("schema_name")

            dself.pop("program")  # not in v1
            dself["protocols"] = self.protocols.convert_v(target_version)

            self_vN = qcmb.models.v1.ManyBodySpecification(**dself)
        else:
            assert False, target_version

        return self_vN


class ManyBodyInput(ProtoModel):
    """Combining the what and how (ManyBodySpecification) with the who (Molecule)."""

    schema_name: Literal["qcschema_many_body_input"] = "qcschema_many_body_input"
    schema_version: Literal[2] = Field(
        2,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    provenance: Provenance = Field(Provenance(**provenance_stamp(__name__)), description=str(Provenance.__doc__))
    id: Optional[str] = None

    specification: ManyBodySpecification = Field(
        ...,
        description="???",
    )
    molecule: Molecule = Field(
        ...,
        description="Target molecule for many-body expansion (MBE) or interaction energy (IE) analysis.",
    )

    def convert_v(
        self, target_version: int, /
    ) -> Union["qcmanybody.models.v1.ManyBodyInput", "qcmanybody.models.v2.ManyBodyInput"]:
        """Convert to instance of particular QCSchema version."""
        import qcmanybody as qcmb

        if check_convertible_version(target_version, error="ManyBodyInput") == "self":
            return self

        dself = self.model_dump()
        if target_version == 1:
            dself.pop("schema_name")
            dself.pop("schema_version")

            dself.pop("id")  # not in v1
            dself.pop("provenance")  # not in v1

            dself["molecule"] = self.molecule.convert_v(target_version)
            dself["specification"] = self.specification.convert_v(target_version)

            self_vN = qcmb.models.v1.ManyBodyInput(**dself)
        else:
            assert False, target_version

        return self_vN


# ====  Properties  =============================================================

# class ManyBodyProperties defined through create_model

manybodyproperties_doc = """
    Named properties of manybody computations following the MolSSI QCSchema.

    All arrays are stored flat but must be reshapable into the dimensions in attribute ``shape``, with abbreviations as follows:

    * nat: number of atoms = :attr:`~qcmanybody.models.v2.ManyBodyProperties.calcinfo_natom`
    * nmc: number of model chemistries = :attr:`~qcmanybody.models.v2.ManyBodyProperties.calcinfo_nmc`
    * nfr: number of fragments = :attr:`~qcmanybody.models.v2.ManyBodyProperties.calcinfo_nfr`
    * nmbe: number of jobs = :attr:`~qcmanybody.models.v2.ManyBodyProperties.calcinfo_nmbe`
    """

MAX_NBODY = int(os.environ.get("QCMANYBODY_MAX_NBODY", 5))  # 5 covers tetramers


json_schema_extras = {
    "energy": {"units": "E_h"},
    "gradient": {"units": "E_h/a0", "shape": ["nat", 3]},
    "Hessian": {"units": "E_h/a0^2", "shape": ["nat" * 3, "nat" * 3]},
}

mbprop = {}

mbprop["schema_name"] = (
    Literal["qcschema_many_body_properties"],
    Field("qcschema_many_body_properties"),
)


# ========  Calcinfo  ===========================================================

mbprop["calcinfo_nmc"] = (
    Optional[int],
    Field(
        None,
        description="The number of model chemistries applied to n-body levels of the computation.",
    ),
)

mbprop["calcinfo_nfr"] = (
    Optional[int],
    Field(
        None,
        description="The number of fragments in the molecule for the computation.",
    ),
)

mbprop["calcinfo_natom"] = (
    Optional[int],
    Field(
        None,
        description="The number of atoms in the computation.",
    ),
)  # alias nat

mbprop["calcinfo_nmbe"] = (
    Optional[int],
    Field(
        None,
        description="The number of real/ghost molecule patterns for the computation.",
    ),
)  # alias NBODY NUMBER

# ========  Canonical  ==========================================================

mbprop["nuclear_repulsion_energy"] = (
    Optional[float],
    Field(
        None,
        description="The nuclear repulsion energy.",
    ),
)

# ret_energy
mbprop["return_energy"] = (
    Optional[float],
    Field(
        None,
        description=f"The interaction energy of the requested method: IE or total (depending on return_total_data) with cp/nocp/vmfc treatment (dep. on first of bsse_type). Always available. Identical to :attr:`~qcelemental.models.ManyBodyResult.return_result` for :attr:`~qcelemental.models.AtomicInput.driver`\\ =\\ :attr:`~qcelemental.models.DriverEnum.energy` computations.",
        json_schema_extra={"units": "E_h"},
    ),
)

# ret_gradient
mbprop["return_gradient"] = (
    Optional[Array[float]],
    Field(
        None,
        description=f"The interaction gradient of the requested method: IE or total (depending on return_total_data) with cp/nocp/vmfc treatment (dep. on first of bsse_type). Available when driver is g/h. Identical to :attr:`~qcelemental.models.ManyBodyResult.return_result` for :attr:`~qcelemental.models.AtomicInput.driver`\\ =\\ :attr:`~qcelemental.models.DriverEnum.gradient` computations.",
        json_schema_extra=json_schema_extras["gradient"],
    ),
)

# ret_hessian
mbprop["return_hessian"] = (
    Optional[Array[float]],
    Field(
        None,
        description=f"The interaction Hessian of the requested method: IE or total (depending on return_total_data) with cp/nocp/vmfc treatment (dep. on first of bsse_type). Available when driver is h. Identical to :attr:`~qcelemental.models.ManyBodyResult.return_result` for :attr:`~qcelemental.models.AtomicInput.driver`\\ =\\ :attr:`~qcelemental.models.DriverEnum.hessian` computations.",
        json_schema_extra=json_schema_extras["Hessian"],
    ),
)

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
            ),
        )

    # CP-CORRECTED TOTAL ENERGY
    mbprop[f"cp_corrected_total_{egh}"] = (
        Optional[typ],
        Field(
            None,
            description=f"Best available total {singular} with cp treatment: cp_corrected_total_{egh}_through_{{max_nbody}}_body. Available when cp in bsse_type & rtd=T{availability_of_derivative}.",
            json_schema_extra=jse,
        ),
    )

    # CP-CORRECTED INTERACTION ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"cp_corrected_interaction_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less 1-body total data for cumulative IE; inputs are total {plural} with cp treatment. Available when when cp in bsse_type & max_nbody>={nb}{availability_of_derivative}. The 1-body quantity is zero by definition.",
                json_schema_extra=jse,
            ),
        )

    # CP-CORRECTED INTERACTION ENERGY
    mbprop[f"cp_corrected_interaction_{egh}"] = (
        Optional[typ],
        Field(
            None,
            description=f"Best available interaction {singular} with cp treatment: cp_corrected_interaction_{egh}_through_{{max_nbody}}_body. Available when cp in bsse_type{availability_of_derivative}.",
            json_schema_extra=jse,
        ),
    )

    # CP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY
    for nb in range(2, MAX_NBODY):
        mbprop[f"cp_corrected_{nb}_body_contribution_to_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less ({nb}-1)-body data for partial IE; inputs are total {plural} w/ cp treat. Available when cp in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ),
        )

    # ========  NOCP E/G/H summary data  ============================================

    # NOCP-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"nocp_corrected_total_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"MBE sum of subsystems of {nb}-body or fewer (cumulative); summed are total {plural} without cp treatment. Available when nocp in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ),
        )

    # NOCP-CORRECTED TOTAL ENERGY
    mbprop[f"nocp_corrected_total_{egh}"] = (
        Optional[typ],
        Field(
            None,
            description=f"Best available total {singular} without cp treatment: nocp_corrected_total_{egh}_through_{{max_nbody}}_body. Available when nocp in bsse_type{availability_of_derivative}.",
            json_schema_extra=jse,
        ),
    )

    # NOCP-CORRECTED INTERATION ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"nocp_corrected_interaction_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less 1-body total data for cumulative IE; inputs are total {plural} without cp treatment. Available when when nocp in bsse_type & max_nbody>={nb}{availability_of_derivative}. The 1-body quantity is zero by definition.",
                json_schema_extra=jse,
            ),
        )

    # TODO note htat TOT 1BODY cp=nocp=vmfc
    # TODO note that summ INTERACTION ENERGY props (w/o explicit -BODY) return 0.0 for max_nbody=1 for completeness

    # NOCP-CORRECTED INTERACTION ENERGY
    mbprop[f"nocp_corrected_interaction_{egh}"] = (
        Optional[typ],
        Field(
            None,
            description=f"Best available interaction {singular} without cp treatment: nocp_corrected_interaction_{egh}_through_{{max_nbody}}_body. Available when nocp in bsse_type{availability_of_derivative}.",
            json_schema_extra=jse,
        ),
    )

    # NOCP-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY
    for nb in range(2, MAX_NBODY):
        mbprop[f"nocp_corrected_{nb}_body_contribution_to_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less ({nb}-1)-body data for partial IE; inputs are total {plural} w/o cp treatment. Available when nocp in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ),
        )

    # ========  VMFC E/G/H summary data  ============================================

    # VMFC-CORRECTED TOTAL ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"vmfc_corrected_total_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"MBE sum of subsystems of {nb}-body or fewer (cumulative); summed are total {plural} with vmfc treatment. Available when vmfc in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ),
        )

    # VMFC-CORRECTED TOTAL ENERGY
    mbprop[f"vmfc_corrected_total_{egh}"] = (
        Optional[typ],
        Field(
            None,
            description=f"Best available total {singular} with vmfc treatment: vmfc_corrected_total_{egh}_through_{{max_nbody}}_body. Available when vmfc in bsse_type{availability_of_derivative}.",
            json_schema_extra=jse,
        ),
    )

    # VMFC-CORRECTED INTERATION ENERGY THROUGH {nb}-BODY
    for nb in range(1, MAX_NBODY):
        mbprop[f"vmfc_corrected_interaction_{egh}_through_{nb}_body"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less 1-body total data for cumulative IE; inputs are total {plural} w/ vmfc treatment. Available when when vmfc in bsse_type & max_nbody>={nb}{availability_of_derivative}. The 1-body quantity is zero by definition.",
                json_schema_extra=jse,
            ),
        )

    # VMFC-CORRECTED INTERACTION ENERGY
    mbprop[f"vmfc_corrected_interaction_{egh}"] = (
        Optional[typ],
        Field(
            None,
            description=f"Best available interaction {singular} with vmfc treatment: vmfc_corrected_interaction_{egh}_through_{{max_nbody}}_body. Available when vmfc in bsse_type{availability_of_derivative}.",
            json_schema_extra=jse,
        ),
    )

    # VMFC-CORRECTED {nb}-BODY CONTRIBUTION TO ENERGY
    for nb in range(2, MAX_NBODY):
        mbprop[f"vmfc_corrected_{nb}_body_contribution_to_{egh}"] = (
            Optional[typ],
            Field(
                None,
                description=f"{nb}-body total data less ({nb}-1)-body total data for partial IE; inputs are total {plural} w/ vmfc treatment. Available when vmfc in bsse_type & max_nbody>={nb}{availability_of_derivative}.",
                json_schema_extra=jse,
            ),
        )


def _validate_arb_max_nbody_fieldnames(cls, values):

    ok_field_name = re.compile(
        # fmt: off
        r"^(?:cp|nocp|vmfc)_corrected_((?:total|interaction)_(?:energy|gradient|hessian)" +
        r"(?:_through_\d+_body)?|\d+_body_contribution_to_(?:energy|gradient|hessian))$"
        # fmt: on
    )

    extra_fields = values.keys() - cls.model_fields.keys()
    baduns = [xtra for xtra in extra_fields if not ok_field_name.match(xtra)]

    if baduns:
        raise ValueError(f"Field names not allowed: {baduns}")

    return values


class ProtoModelSkipDefaults(ProtoModel):

    # fields filtered in model_validator
    model_config = ExtendedConfigDict(serialize_skip_defaults=True, force_skip_defaults=True, extra="allow")


if TYPE_CHECKING:
    ManyBodyProperties = ProtoModelSkipDefaults
else:
    # if/else suppresses a warning about using a dynamically generated class as Field type in ManyBodyResults
    # * deprecated but works: root_validator(skip_on_failure=True)(_validate_arb_max_nbody_fieldnames)
    ManyBodyProperties = create_model(
        "ManyBodyProperties",
        # __doc__=manybodyproperties_doc,  # needs later pydantic
        __base__=ProtoModelSkipDefaults,
        __validators__={"validator1": model_validator(mode="before")(_validate_arb_max_nbody_fieldnames)},
        **mbprop,
    )


def _qcvars_translator(cls, reverse: bool = False) -> Dict[str, str]:
    """Form translation map between many-body results QCSchema and Psi4/QCDB terminologies.

    Parameters
    ----------
    reverse
        Keys are QCVariable names (`reverse=True`) rather than QCSchema names (default; `reverse=False`).

    Returns
    -------
    dict
        Map from ManyBodyProperties field names to QCVariable names, or reverse.

    """
    qcvars_to_mbprop = {}
    for skprop in ManyBodyProperties.model_fields.keys():
        qcvar = skprop.replace("_body", "-body").replace("_corr", "-corr").replace("_", " ").upper()
        qcvars_to_mbprop[qcvar] = skprop
    for ret in ["energy", "gradient", "hessian"]:
        qcvars_to_mbprop[f"CURRENT {ret.upper()}"] = f"return_{ret}"
    qcvars_to_mbprop["NBODY NUMBER"] = "calcinfo_nmbe"

    if reverse:
        return qcvars_to_mbprop
    else:
        return {v: k for k, v in qcvars_to_mbprop.items()}


ManyBodyProperties.to_qcvariables = classmethod(_qcvars_translator)


# ====  Results  ================================================================


class ManyBodyResult(SuccessfulResultBase):

    schema_name: Literal["qcschema_many_body_result"] = "qcschema_many_body_result"
    schema_version: Literal[2] = Field(
        2,
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
    success: Literal[True] = Field(
        True,
        description="A boolean indicator that the operation succeeded or failed. Allows programmatic assessment of "
        "all results regardless of if they failed or succeeded by checking `result.success`.",
    )

    # native_files placeholder for when any mbe programs supply extra files or need an input file. no protocol at present
    native_files: Dict[str, Any] = Field({}, description="DSL files.")

    molecule: Molecule = Field(
        ...,
        description="The molecule in results fragmentation and frame. Since QCManyBody doesn't disrupt the mol, should be identical to input_data.molecule.",
    )

    properties: ManyBodyProperties = Field(..., description=str(ManyBodyProperties.__doc__))
    cluster_properties: Dict[str, AtomicProperties] = Field(
        ...,
        description="The key results for each subsystem species computed. Keys contain modelchem, real and ghost information (e.g., `'[\"(auto)\", [2], [1, 2, 3]]'`). Values are total e/g/H/property results. Array values, if present, are sized and shaped for the full supersystem.",
    )
    cluster_results: Dict[str, AtomicResult] = Field(..., description="Detailed results")
    return_result: Union[float, Array[float], Dict[str, Any]] = Field(
        ...,
        description="The primary return specified by the :attr:`~qcelemental.models.AtomicInput.driver` field. Scalar if energy; array if gradient or hessian; dictionary with property keys if properties.",
    )
    stdout: Optional[str] = Field(
        None,
        description="The primary logging output of the program, whether natively standard output or a file. Presence vs. absence (or null-ness?) configurable by protocol.",
    )
    stderr: Optional[str] = Field(None, description="The standard error of the program execution.")

    @field_validator("cluster_results")
    def _cluster_results(cls, value, values: ValidationInfo):
        crp = values.data["input_data"].specification.protocols.cluster_results
        if crp == "all":
            return value
        elif crp == "none":
            return {}
        else:
            raise ValueError(f"Protocol `component_resutls:{crp}` is not understood")

    def convert_v(
        self, target_version: int, /
    ) -> Union["qcmanybody.models.v1.ManyBodyResult", "qcmanybody.models.v2.ManyBodyResult"]:
        """Convert to instance of particular QCSchema version."""
        import qcmanybody as qcmb

        if check_convertible_version(target_version, error="ManyBodyResult") == "self":
            return self

        dself = self.model_dump()
        if target_version == 1:
            dself.pop("schema_name")  # changed in v1
            dself.pop("schema_version")  # changed in v1

            # for input_data, work from model, not dict, to use convert_v
            dself["input_data"] = self.input_data.convert_v(1).model_dump()  # exclude_unset=True, exclude_none=True

            dself.pop("native_files")
            dself.pop("molecule")

            dself["component_properties"] = dself.pop("cluster_properties")
            dself["component_results"] = {
                k: atres.convert_v(target_version) for k, atres in self.cluster_results.items()
            }
            dself.pop("cluster_results")

            self_vN = qcmb.models.v1.ManyBodyResult(**dself)
        else:
            assert False, target_version

        return self_vN
