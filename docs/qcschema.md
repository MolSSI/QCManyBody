
<!-- ====  Inputs  ================================================================= -->

::: qcmanybody.models.BsseEnum
    options:
        show_root_heading: true
        show_source: true


::: qcmanybody.models.ManyBodyKeywords
    options:
        show_root_heading: true

$pydantic: qcmanybody.models.manybody_input_pydv1.ManyBodyKeywords


::: qcmanybody.models.ManyBodySpecification
    options:
        show_root_heading: true

$pydantic: qcmanybody.models.manybody_input_pydv1.ManyBodySpecification


::: qcmanybody.models.ManyBodyInput
    options:
        show_root_heading: true

$pydantic: qcmanybody.models.manybody_input_pydv1.ManyBodyInput


<!-- ====  Protocols  ============================================================== -->

<!--
::: qcmanybody.models.manybody_pydv1.ManyBodyProtocolEnum


::: qcmanybody.models.manybody_input_pydv1.ManyBodyProtocols
    options:
        show_root_heading: true

$pydantic: qcmanybody.models.manybody_input_pydv1.ManyBodyProtocols
-->


<!-- ====  Properties/Outputs  ===================================================== -->

!!! note

    The properties model is generated dynamically based on a constant
    ``MAX_NBODY``. To not overload the docs table, this is set to 5, which
    covers full calculations on tetramers. In practice this isn't a problem for
    larger clusters because ``cp_corrected_total_energy_through_12_body``, for
    example, is allowed dynamically for a model instance. Nevertheless, to use a
    larger ``ManyBodyKeywords.max_nbody``, reset this value *outside* the interpreter.

        python -c "import qcmanybody as qcmb;print(qcmb.models.MAX_NBODY)"
        #> 5
        export QCMANYBODY_MAX_NBODY=9  # explicitly enumerates octamer properties
        python -c "import qcmanybody as qcmb;print(qcmb.models.MAX_NBODY)"
        #> 9


::: qcmanybody.models.ManyBodyResultProperties
    options:
        show_root_heading: true

$pydantic: qcmanybody.models.ManyBodyResultProperties


::: qcmanybody.models.ManyBodyResult
    options:
        show_root_heading: true

$pydantic: qcmanybody.models.ManyBodyResult


<!-- ====  Misc.  ================================================================== -->

<!-- $pydantic: qcmanybody.models.manybody_pydv1.AtomicSpecification -->
<!--
AtomicSpecification
ResultsBase
SuccessfulResultBase
-->

<!--
    options:
        merge_init_into_class: false
        group_by_category: false
        # explicit members list so we can set order and include `__init__` easily
        members:
          - __init__
          - molecule
          - model_config
          - model_computed_fields
          - model_extra
          - model_fields
          - model_fields_set
          - model_construct
          - model_copy
          - model_dump
          - model_dump_json
          - model_json_schema
          - model_parametrized_name
          - model_post_init
          - model_rebuild
          - model_validate
          - model_validate_json
          - copy
-->
