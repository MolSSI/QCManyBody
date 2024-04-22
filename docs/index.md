# QCManyBody Documentation

QCManyBody is a python package for running quantum chemistry manybody expansions and interaction calculations in a
package-independent way.

## Installation

Currently, the package is only available on [GitHub](https://github.com/MolSSI/QCManyBody). To install directly
from GitHub, use the following command:

```bash
pip install git+https://github.com/MolSSI/QCManyBody.git
```

## Package Overview

The package has two main interfaces. The high-level interface allows for a comprehensive workflow, where the user
provides complete information about the calculation, including the full specification (method, basis set, etc.) of the
manybody calculation. This is designed to work with [QCEngine](https://github.com/MolSSI/QCEngine) or other packages
that implement the [QCSchema](https://github.com/MolSSI/QCSchema).

For more information, see (High-level interface)(high-level-interface.md).

QCManyBody also has a core low-level interface that allows for more flexibility in how the calculations are run. This
interface generally takes a molecule and an arbitrary definition of quantum chemistry specifications, and expects
the user to run them themselves.

For more information, see (Core interface)(core-interface.md).

## Table of Contents

1. [Common Keywords and Options](keywords.md)
2. [High-level Interface](high-level-interface.md)
3. [Core Interface](core-interface.md)
4. [Results](results.md)
5. [How-To Guides](how-to-guides.md)