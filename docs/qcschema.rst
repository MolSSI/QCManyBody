========
QCSchema
========

Inputs
======

.. automodule:: qcmanybody.models.v1
   :members: BsseEnum, ManyBodyKeywords, ManyBodyProtocols, ManyBodySpecification, ManyBodyInput
   :undoc-members:
   :show-inheritance:


Properties/Outputs
==================

.. note::

    The properties model is generated dynamically based on a constant
    ``MAX_NBODY``. To not overload the docs table, this is set to 5, which
    covers full calculations on tetramers. In practice this isn't a problem for
    larger clusters because ``cp_corrected_total_energy_through_12_body``, for
    example, is allowed dynamically for a model instance. Nevertheless, to use a
    larger ``ManyBodyKeywords.max_nbody``, reset this value *outside* the interpreter.

    .. code-block:: python

        python -c "import qcmanybody as qcmb;print(qcmb.models.MAX_NBODY)"
        #> 5
        export QCMANYBODY_MAX_NBODY=9  # explicitly enumerates octamer properties
        python -c "import qcmanybody as qcmb;print(qcmb.models.MAX_NBODY)"
        #> 9


.. automodule:: qcmanybody.models.v1
   :members: ManyBodyResultProperties, ManyBodyResult
   :undoc-members:
   :show-inheritance:
