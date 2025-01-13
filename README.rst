.. _images:

.. image:: doc/source/_static/parsnip_header_dark.svg
  :width: 600


.. _header:

..
  TODO: set up Readthedocs, PyPI, and conda-forge

|ReadTheDocs|
|PyPI|
|conda-forge|

.. |ReadTheDocs| image:: https://readthedocs.org/projects/parsnip-cif/badge/?version=latest
   :target: http://parsnip-cif.readthedocs.io/en/latest/?badge=latest
.. |PyPI| image:: https://img.shields.io/pypi/v/parsnip-cif.svg
   :target: https://pypi.org/project/parsnip-cif/
.. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/parsnip-cif.svg
   :target: https://anaconda.org/conda-forge/parsnip-cif


.. _introduction:

**parsnip** is a minimal Python library for parsing `CIF <https://www.iucr.org/resources/cif>`_ files. While its primary focus is on simplicity and portability, performance-oriented design choices are made where possible.

.. _parse:

Importing ``parsnip`` allows users to read `CIF 1.1 <https://www.iucr.org/resources/cif/spec/version1.1>`_ files, as well as many features from the `CIF 2.0 <https://www.iucr.org/resources/cif/cif2>`_ and `mmCIF <https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner’s-guide-to-pdb-structures-and-the-pdbx-mmcif-format>`_ formats.
Creating a :class:`~.CifFile` object provides easy access to name-value :attr:`~.CifFile.pairs`, as well
as `loop\_`-delimited :attr:`~.CifFile.tables`. Data entries can be extracted as python primitives or
numpy arrays for further use.

.. _installing:

Setup
-----

**parsnip** may be installed with **pip** or from **conda-forge**.


Installation via pip
^^^^^^^^^^^^^^^^^^^^

.. code:: bash

    python -m pip install parsnip

Installation via conda-forge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash

    conda install -c conda-forge parsnip


Installation from source
^^^^^^^^^^^^^^^^^^^^^^^^

First, clone the repository:

.. code:: bash

    git clone https://github.com/glotzerlab/parsnip.git
    cd parsnip

Then, choose one of the following. While **parsnip** is only dependent on Numpy,
additional dependencies are required to run the tests and build the docs.

.. code:: bash

    pip install .            # Install with no additional dependencies
    pip install .[tests]     # Install with dependencies required to run tests
    pip install .[tests,doc] # Install with dependencies required to run tests and make docs

Dependencies
^^^^^^^^^^^^

.. code:: text

   numpy>=1.19
   more-itertools

.. _contributing:
