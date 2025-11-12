# QCManyBody Documentation

QCManyBody is a python package for running quantum chemistry many-body expansions and interaction calculations in a
package-independent way.

## Installation

QCManyBody is available from [PyPI](https://pypi.org/project/qcmanybody) and from
[conda-forge](https://anaconda.org/conda-forge/qcmanybody).

```bash
# Installation from PyPI
pip install qcmanybody

# Installation from conda-forge
conda install -c conda-forge qcmanybody
```

To install the latest development version directly from
[GitHub](https://github.com/MolSSI/QCManyBody), you can use the following command:

```bash
pip install git+https://github.com/MolSSI/QCManyBody.git
```

## Package Overview

The package has two main interfaces. The high-level interface allows for a comprehensive workflow, where the user
provides complete information about the calculation, including the full specification (method, basis set, etc.) of the
many-body calculation. This is designed to work with [QCEngine](https://github.com/MolSSI/QCEngine) or other packages
that implement the [QCSchema](https://github.com/MolSSI/QCSchema).

For more information, see [High-level interface](high-level-interface.md).

QCManyBody also has a low-level _core_ interface that allows for more flexibility in how the calculations are run. This
interface generally takes a molecule and an arbitrary definition of quantum chemistry specifications, and expects
the user to run them themselves.

For more information, see [Core interface](core-interface.md).

QCManyBody supports **parallel execution** for many-body calculations, enabling significant speedups by running independent
quantum chemistry tasks concurrently. The parallel execution module provides multiple backends including multiprocessing
for single-node parallelism and MPI for distributed HPC clusters.

For more information, see [Parallel Execution Guide](parallel_execution_guide.md).

## Table of Contents

1. [Common Keywords and Options](keywords.md)
2. [High-level Interface](high-level-interface.md)
3. [Core Interface](core-interface.md)
4. [Command-Line Interface](cli_guide.md)
5. [Results](results.md)
6. [How-To Guides](how-to-guides.md)

### Parallel Execution

7. [Parallel Execution Guide](parallel_execution_guide.md) - Getting started with parallel calculations
8. [Parallel API Reference](parallel_api_reference.md) - Complete API documentation
9. [Parallel Migration Guide](parallel_migration_guide.md) - Migrating existing code to use parallel execution

## [Changelog](changelog.md)

<!--

PYTHONPATH=docs/extensions:. mkdocs serve

python -m sphinx.ext.intersphinx 'https://molssi.github.io/QCManyBody/objects.inv'
python -m sphinx.ext.intersphinx 'http://127.0.0.1:8000/QCManyBody/objects.inv'

mkdocs build

-->
