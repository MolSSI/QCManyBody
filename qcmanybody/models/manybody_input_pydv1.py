from warnings import warn

import qcmanybody

_qcsk_v2_default_v1_importpathschange = "0.70.0"  # qcelemental.models.common_models

_nonapi_file = "manybody_input_pydv1"

warn(
    f"qcmanybody.models.{_nonapi_file} should be accessed through qcmanybody.models (or qcmanybody.models.v1 or .v2 for fixed QCSchema version). The 'models.{_nonapi_file}' route will be removed as soon as v{_qcsk_v2_default_v1_importpathschange}.",
    DeprecationWarning,
)

BsseEnum = qcmanybody.models.v1.BsseEnum
ManyBodyProtocols = qcmanybody.models.v1.ManyBodyProtocols
ManyBodyKeywords = qcmanybody.models.v1.ManyBodyKeywords
ManyBodySpecification = qcmanybody.models.v1.ManyBodySpecification
ManyBodyInput = qcmanybody.models.v1.ManyBodyInput
