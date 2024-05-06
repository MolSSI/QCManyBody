from typing import Annotated, List, Literal, Union

try:
    from pydantic.v1 import Field
except ImportError:
    from pydantic import Field

from qcelemental.models import AtomicResult, OptimizationInput, OptimizationResult
from qcelemental.models.procedures import QCInputSpecification

from .manybody_input_pydv1 import ManyBodySpecification
from .manybody_output_pydv1 import ManyBodyResult


# note that qcel AtomicResult.schema_name needs editing
ResultTrajectories = Annotated[Union[AtomicResult, ManyBodyResult], Field(discriminator='schema_name')]

class GeneralizedOptimizationInput(OptimizationInput):
    schema_name: Literal["qcschema_generalizedoptimizationinput"] = "qcschema_generalizedoptimizationinput"
    schema_version: int = 1
    input_specification: Union[QCInputSpecification, ManyBodySpecification] = Field(..., description="ordinary or mbe grad spec")


class GeneralizedOptimizationResult(OptimizationResult):
    schema_name: Literal["qcschema_generalizedoptimizationresult"] = "qcschema_generalizedoptimizationresult"
    trajectory: List[ResultTrajectories] = Field(
        ..., description="A list of ordered Result objects for each step in the optimization."
    )
    input_specification: Union[QCInputSpecification, ManyBodySpecification] = Field(..., description="ordinary or mbe grad spec")
