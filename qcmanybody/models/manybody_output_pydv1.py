from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING, Any, Dict, Literal, Optional, Union

# v2: from pydantic import create_model, Field, field_validator, FieldValidationInfo
try:
    from pydantic.v1 import Field, create_model, root_validator, validator
except ImportError:
    from pydantic import create_model, Field, validator, root_validator

from qcelemental.models import ProtoModel, Provenance
from qcelemental.models.results import AtomicResult, AtomicResultProperties
from qcelemental.models.types import Array

from .manybody_input_pydv1 import ManyBodyInput, SuccessfulResultBase

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

MAX_NBODY = int(os.environ.get("QCMANYBODY_MAX_NBODY", 5))  # 5 covers tetramers


json_schema_extras = {
    "energy": {"units": "E_h"},
    "gradient": {"units": "E_h/a0", "shape": ["nat", 3]},
    "Hessian": {"units": "E_h/a0^2", "shape": ["nat" * 3, "nat" * 3]},
}

mbprop = {}

mbprop["schema_name"] = (
    Literal["qcschema_manybodyproperties"],
    Field("qcschema_manybodyproperties"),
)
mbprop["schema_version"] = (
    Literal[1],
    Field(
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    ),
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

    extra_fields = values.keys() - cls.__fields__.keys()
    baduns = [xtra for xtra in extra_fields if not ok_field_name.match(xtra)]

    if baduns:
        raise ValueError(f"Field names not allowed: {baduns}")

    return values


class ProtoModelSkipDefaults(ProtoModel):

    class Config(ProtoModel.Config):
        serialize_skip_defaults = True
        force_skip_defaults = True
        extra: str = "allow"  # fields filtered in root_validator


if TYPE_CHECKING:
    ManyBodyResultProperties = ProtoModelSkipDefaults
else:
    # if/else suppresses a warning about using a dynamically generated class as Field type in ManyBodyResults
    ManyBodyResultProperties = create_model(
        "ManyBodyResultProperties",
        # __doc__=manybodyresultproperties_doc,  # needs later pydantic
        __base__=ProtoModelSkipDefaults,
        __validators__={"validator1": root_validator(_validate_arb_max_nbody_fieldnames)},
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
        Map from ManyBodyResultProperties field names to QCVariable names, or reverse.

    """
    qcvars_to_mbprop = {}
    # v2: for skprop in ManyBodyResultProperties.model_fields.keys():
    for skprop in cls.__fields__.keys():
        qcvar = skprop.replace("_body", "-body").replace("_corr", "-corr").replace("_", " ").upper()
        qcvars_to_mbprop[qcvar] = skprop
    for ret in ["energy", "gradient", "hessian"]:
        qcvars_to_mbprop[f"CURRENT {ret.upper()}"] = f"return_{ret}"
    qcvars_to_mbprop["NBODY NUMBER"] = "calcinfo_nmbe"

    if reverse:
        return qcvars_to_mbprop
    else:
        return {v: k for k, v in qcvars_to_mbprop.items()}


ManyBodyResultProperties.to_qcvariables = classmethod(_qcvars_translator)


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
    component_properties: Dict[str, AtomicResultProperties] = Field(
        ...,
        description="The key results for each subsystem species computed. Keys contain modelchem, real and ghost information (e.g., `'[\"(auto)\", [2], [1, 2, 3]]'`). Values are total e/g/H/property results. Array values, if present, are sized and shaped for the full supersystem.",
    )
    component_results: Dict[str, AtomicResult] = Field({}, description="Detailed results")
    return_result: Union[float, Array[float], Dict[str, Any]] = Field(
        ...,
        description="The primary return specified by the :attr:`~qcelemental.models.AtomicInput.driver` field. Scalar if energy; array if gradient or hessian; dictionary with property keys if properties.",
    )
    stdout: Optional[str] = Field(
        None,
        description="The primary logging output of the program, whether natively standard output or a file. Presence vs. absence (or null-ness?) configurable by protocol.",
    )
    stderr: Optional[str] = Field(None, description="The standard error of the program execution.")
    # v2: success: Literal[True] = Field(True, description="Always `True` for a successful result")

    @validator("component_results")
    def _component_results(cls, value, values):
        crp = values["input_data"].specification.protocols.component_results
        if crp == "all":
            return value
        elif crp == "none":
            return {}
        else:
            raise ValueError(f"Protocol `component_resutls:{crp}` is not understood")
