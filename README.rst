.. _header:

.. image:: _static/parsnip_header_dark.svg
  :width: 600
  :class: only-light

.. image:: _static/parsnip_header_light.svg
  :width: 600
  :class: only-dark

..
  TODO: set up Readthedocs, PyPI, and conda-forge

|ReadTheDocs|
|PyPI|
|conda-forge|

.. |ReadTheDocs| image:: https://readthedocs.org/projects/parsnip/badge/?version=latest
   :target: http://parsnip.readthedocs.io/en/latest/?badge=latest
.. |PyPI| image:: https://img.shields.io/pypi/v/parsnip.svg
   :target: https://pypi.org/project/parsnip/
.. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/parsnip.svg
   :target: https://anaconda.org/conda-forge/parsnip


**parsnip** is a minimal Python library for parsing `CIF <https://www.iucr.org/resources/cif>`_ files. While its primary focus is on simplicity and portability, performance-oriented design choices are made where possible.


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

Install Requirements
^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../requirements.txt
  :language: text
