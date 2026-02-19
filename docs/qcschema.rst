========
QCSchema
========

.. note:: QCManyBody imports *without* QCSchema/Pydantic version specification (
   ``qcmanybody.models`` vs ``models.v1`` and ``models.v2`` and ``qcmanybody.computer`` vs. ``v1.computer`` and ``v2.computer)
   are *always* pointing to v1 for continuity. This will be true until QCManyBody depends on QCElemental v0.70.0.
   For technical constraints, these docs always reference v2 forms. See :docs:`v0.5.2` for QCSchema v1-based documentation.

Inputs
======

.. automodule:: qcmanybody.models.v2
   :members: BsseEnum, ManyBodyKeywords, ManyBodyProtocols, ManyBodySpecification, ManyBodyInput
   :undoc-members:
   :show-inheritance:
   :noindex:


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

        python -c "import qcmanybody as qcmb;print(qcmb.models.v2.MAX_NBODY)"
        #> 5
        export QCMANYBODY_MAX_NBODY=9  # explicitly enumerates octamer properties
        python -c "import qcmanybody as qcmb;print(qcmb.models.v2.MAX_NBODY)"
        #> 9


.. automodule:: qcmanybody.models.v2
   :members: ManyBodyProperties, ManyBodyResult
   :undoc-members:
   :show-inheritance:
   :noindex:
