=========================
QCManyBody Documentation
=========================

QCManyBody is a python package for running quantum chemistry many-body expansions and interaction calculations in a
package-independent way.

Installation
============

QCManyBody is available from `PyPI <https://pypi.org/project/qcmanybody>`_ and from
`conda-forge <https://anaconda.org/conda-forge/qcmanybody>`_.

.. code-block:: bash

    # Installation from PyPI
    pip install qcmanybody

    # Installation from conda-forge
    conda install -c conda-forge qcmanybody

To install the latest development version directly from
`GitHub <https://github.com/MolSSI/QCManyBody>`_, you can use the following command:

.. code-block:: bash

    pip install git+https://github.com/MolSSI/QCManyBody.git

Package Overview
================

The package has two main interfaces. The high-level interface allows for a comprehensive workflow, where the user
provides complete information about the calculation, including the full specification (method, basis set, etc.) of the
many-body calculation. This is designed to work with `QCEngine <https://github.com/MolSSI/QCEngine>`_ or other packages
that implement the `QCSchema <https://github.com/MolSSI/QCSchema>`_.

For more information, see :doc:`high-level-interface`.

QCManyBody also has a low-level *core* interface that allows for more flexibility in how the calculations are run. This
interface generally takes a molecule and an arbitrary definition of quantum chemistry specifications, and expects
the user to run them themselves.

For more information, see :doc:`core-interface`.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   high-level-interface
   core-interface
   qcschema
   api
   keywords
   results
   changelog

.. python -m sphinx.ext.intersphinx 'https://molssi.github.io/QCManyBody/objects.inv'
.. python -m sphinx.ext.intersphinx 'http://127.0.0.1:8000/QCManyBody/objects.inv'
