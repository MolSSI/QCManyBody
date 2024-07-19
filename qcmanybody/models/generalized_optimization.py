from typing import Annotated, List, Literal, Union

try:
    from pydantic.v1 import Field
except ImportError:
    from pydantic import Field

from qcelemental.models import AtomicResult, OptimizationInput, OptimizationResult
from qcelemental.models.procedures import QCInputSpecification

from .manybody_input_pydv1 import ManyBodySpecification
from .manybody_output_pydv1 import ManyBodyResult

# note that qcel QCInputSpecification and AtomicResult.schema_name needs editing
ResultTrajectories = Annotated[
    Union[AtomicResult, ManyBodyResult],
    Field(
        discriminator="schema_name",
        description="A result object for a single step in the optimization. Either an ordinary atomic/single-point or a many-body result.",
    ),
]
InputSpecifications = Annotated[
    Union[QCInputSpecification, ManyBodySpecification],
    Field(
        discriminator="schema_name",
        description="A directive to compute a gradient. Either an ordinary atomic/single-point or a many-body spec.",
    ),
]


class GeneralizedOptimizationInput(OptimizationInput):
    schema_name: Literal["qcschema_generalizedoptimizationinput"] = "qcschema_generalizedoptimizationinput"
    schema_version: int = 1
    input_specification: InputSpecifications


class GeneralizedOptimizationResult(OptimizationResult):
    schema_name: Literal["qcschema_generalizedoptimizationresult"] = "qcschema_generalizedoptimizationresult"
    trajectory: List[ResultTrajectories] = Field(
        ..., description="A list of ordered Result objects for each step in the optimization."
    )
    input_specification: InputSpecifications
