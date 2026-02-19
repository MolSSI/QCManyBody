=====================
High-Level Interface
=====================

The high-level interface takes advantage of :doc:`qcschema`
data structures for input and output and defines an easy way to run
the QC computations through QCEngine.
The high-level interface uses the core-interface under the hood.
See the :class:`~qcmanybody.models.v2.ManyBodyInput` QCSchema class for the input options :class:`~qcmanybody.models.v2.ManyBodyKeywords`.
Details on forming fragmented ``Molecule`` s are at
:ref:`molecule input <sec:keywords_molecule>`

The high-level interface defines a :class:`~qcmanybody.v2.computer.ManyBodyComputer` class in
computer.py whose input is :class:`~qcmanybody.models.v2.ManyBodyInput` schema and whose output is
:class:`~qcmanybody.models.v2.ManyBodyResult` schema.
It provides :func:`qcmanybody.v2.computer.ManyBodyComputer.from_manybodyinput()` to run through
QCEngine.
See ``test_highlevel_interface_example()`` in
`test_example.py <https://github.com/MolSSI/QCManyBody/blob/main/qcmanybody/tests/test_examples.py>`_
for a working example.
The :class:`~qcmanybody.models.v2.ManyBodyComputer` provides the missing link from the core interface
to run QC computations. A strategy that Psi4 uses is to define its own
``ManyBodyComputer`` (inheriting from the high-level interface) that
runs QC computations the way Psi4 wants to while still using all the
validation and structure provided by the QCSchema I/O classes.

.. note:: QCManyBody imports *without* QCSchema/Pydantic version specification (
   models vs ``models.v1`` and ``models.v2`` and ``computer`` vs. ``v1.computer`` and ``v2.computer)
   are *always* pointing to v1 for continuity. This will be true until QCManyBody depends on QCElemental v0.70.0.
   For technical constraints, these docs always reference v2 forms. See :docs:`v0.5.2` for QCSchema v1-based documentation.
