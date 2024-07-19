# QCManyBody Testing

Run tests with pytest or see [GHA](/.github/workflows/ci.yml) for an example. Individual tests will skip if they need
requirements (QCEngine, QC runner, zstandard, etc.) not detected in the current environment.

```bash
pytest -v qcmanybody/
```

There are inadvertantly a few classes of tests:

### Examples

* includes [test_examples](test_examples.py)
* tests core and high-level interfaces
* easy to read
* as set up, one test runs w/o any QC runners, while the other two tests require Psi4 and NWChem.
* _A_ way to run the examples is to create a conda environment from the file below (Linux/Mac/Win), activate it, then:
  * run it through pytest as `pytest -v qcmanybody/ -k "examples"`
  * -or- run it through python as `python qcmanybody/tests/test_examples.py`

```yaml
name: test
channels:
  - conda-forge
  - nodefaults
dependencies:
  - python
  - qcmanybody
  - qcengine
  - pytest
  - psi4=1.9.1
  - nwchem
```

### Unit Tests

* includes [test_schema_keywords](test_schema_keywords.py), [test_utils](test_utils.py), [test_core_auxiliary](test_core_auxiliary.py)
* runs quickly with no extra requirements
* tests elements of QCManyBody like input schema and utilities

### Static-Data Regression Tests

* includes [test_core_singlelevel](test_core_singlelevel.py), [test_core_multilevel](test_core_multilevel.py), [test_core_multilevel_ss](test_core_multilevel_ss.py)
* tests core interface
* HARD TO READ! These test many combinations of inputs and a dictionary comparison against pre-stored outputs, so parts
  are scattered. Users looking for examples aren't advised to look at these.
* requires zstandard to be installed
* no QC runners needed. Output data from Psi4 are stored and may need to be regenerated as the interface changes.
* no QC calculations are run, so testing is instantaneous.

### End-to-End Correctness Tests

* includes [test_computer_he4_singlelevel](test_computer_he4_singlelevel.py), [test_computer_he4_multilevel](test_computer_he4_multilevel.py), [test_computer_het4_gradient.py](test_computer_het4_gradient.py)
* tests high-level (and through it, core) interfaces
* HARD TO READ! These test many combinations of inputs and many details of outputs, so to avoid repetition, there is
  much parameterization and dictionary lookups. Users looking for examples aren't advised to look at these.
* requires QCEngine be installed
* requires one of Psi4, NWChem, or Cfour be installed (actual QCManyBody tests are the same with each QC runner)
* small QC calculations are actually run, so takes a few minutes. Prune with `-k "not (3b or 4b)"`, for example.
