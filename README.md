QCManyBody
==========

<p align="left">
    <picture>
    <img alt="QCManyBody Logo" src="https://github.com/MolSSI/QCManyBody/blob/main/docs/logo.png" height="150px">
    </picture>
</p>

[![GitHub Actions](https://img.shields.io/github/actions/workflow/status/MolSSI/QCManyBody/ci.yml?logo=github)](https://github.com/MolSSI/QCManyBody/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/loriab/QCManyBody/graph/badge.svg?token=E4S0706HJ0)](https://codecov.io/gh/loriab/QCManyBody)
[![Documentation Status](https://img.shields.io/github/actions/workflow/status/MolSSI/QCManyBody/ci.yml?label=docs&logo=readthedocs&logoColor=white)](https://molssi.github.io/QCManyBody/)
[![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/qcmanybody?color=blue&logo=anaconda&logoColor=white)](https://anaconda.org/conda-forge/qcmanybody)
![python](https://img.shields.io/badge/python-3.8+-blue.svg)

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

## Documentation

Full documentation is available at [https://molssi.github.io/QCManyBody/](https://molssi.github.io/QCManyBody/)

## Authors

* Asem Alenaizan, [@alenaizan](https://github.com/alenaizan), original Psi4 implementations of vmfc Hessians, multi-level, and embedded point charges
* Lori A. Burns, [@loriab](https://github.com/loriab), ManyBody QCSchema and high-level interface
* Benjamin P. Pritchard, [@bennybp](https://github.com/bennybp), core interface and QCArchive integration
* Daniel G. A. Smith, [@dgasmith](https://github.com/dgasmith), original Psi4 implementations of nocp, cp, and vmfc single-level e/g/H and distributed driver integration

## Citation [![doi](https://img.shields.io/badge/doi-10.1063/5.0231843-5077AB.svg)](https://doi.org/10.1063/5.0231843)

The journal article reference describing QCManyBody is:

L. A. Burns, C. D. Sherrill, B. P. Pritchard,
"QCManyBody: A Flexible Implementation of the Many-Body Expansion",
J. Chem. Phys. 161(15) 152501 (2024).

## Demonstration

An example of the core and high-level interfaces can be found in [test_examples](qcmanybody/tests/test_examples.py) with
directions in [tests/README](qcmanybody/tests/README.md).
