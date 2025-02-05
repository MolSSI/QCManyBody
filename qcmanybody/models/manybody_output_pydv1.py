from warnings import warn

from qcelemental.models.common_models import _qcsk_v2_default_v1_importpathschange

_nonapi_file = "manybody_output_pydv1"

warn(
    f"qcmanybody.models.{_nonapi_file} should be accessed through qcmanybody.models (or qcmanybody.models.v1 or .v2 for fixed QCSchema version). The 'models.{_nonapi_file}' route will be removed as soon as v{_qcsk_v2_default_v1_importpathschange}.",
    DeprecationWarning,
)

ManyBodyResultProperties = qcmanybody.models.v1.ManyBodyResultProperties
ManyBodyResult = qcmanybody.models.v1.ManyBodyResult
