::: qcmanybody.ManyBodyCore
    options:
        show_root_heading: true

::: qcmanybody.ManyBodyComputer
    options:
        show_root_heading: true
        inherited_members: true
        members:
          - from_manybodyinput

$pydantic: qcmanybody.computer.ManyBodyComputer

::: qcmanybody.utils
    options:
        show_root_heading: true

::: qcmanybody.builder
    options:
        show_root_heading: true
