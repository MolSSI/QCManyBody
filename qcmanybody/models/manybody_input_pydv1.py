from warnings import warn

from qcelemental.models.common_models import _qcsk_v2_default_v1_importpathschange

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
