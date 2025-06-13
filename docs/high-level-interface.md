# High-Level Interface

The high-level interface takes advantage of [QCSchema](qcschema.md)
data structures for input and output and defines an easy way to run
the QC computations through QCEngine.
The high-level interface uses the core-interface under the hood.
See the `ManyBodyInput` QCSchema class for the input options `ManyBodyKeywords.
Details on forming fragmented `Molecule`s are at
[moleule input](keywords.md#molecule)

The high-level interface defines a `ManyBodyComputer` class in
computer.py whose input is `ManyBodyInput` schema and whose output is
`ManyBodyResult` schema.
It provides `ManyBodyComputer.from_manybodyinput()` to run through
QCEngine.
See `test_highlevel_interface_example()` in
[test_example.py](https://github.com/MolSSI/QCManyBody/blob/main/qcmanybody/tests/test_examples.py)
for a working example.
The `ManyBodyComputer` provides the missing link from the core interface
to run QC computations. A strategy that Psi4 uses is to define its own
`ManyBodyComputer` (inheriting from the high-level interface) that
runs QC computations the way Psi4 wants to while still using all the
validation and structure provided by the QCSchema I/O classes.
