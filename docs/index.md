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

## Table of Contents

1. [Common Keywords and Options](keywords.md)
2. [High-level Interface](high-level-interface.md)
3. [Core Interface](core-interface.md)
4. [Results](results.md)
5. [How-To Guides](how-to-guides.md)

<!--

PYTHONPATH=docs/extensions:. mkdocs serve

python -m sphinx.ext.intersphinx 'https://molssi.github.io/QCManyBody/objects.inv'
python -m sphinx.ext.intersphinx 'http://127.0.0.1:8000/QCManyBody/objects.inv'

mkdocs build

-->
