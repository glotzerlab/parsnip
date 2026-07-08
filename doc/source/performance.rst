.. _performance:

===========
Performance
===========

Although ``parsnip`` is optimized for accuracy rather than performace, our parsing
strategy results in file reads and unit-cell reconstructions that are 3-500x faster than
comparable tools. The following data makes use of
`ltalirz/cif-parsing-benchmark <https://github.com/ltalirz/cif-parsing-benchmark>`_,
a benchmark that profiles the parsing throughput of several open-source CIF libraries.

.. figure:: ../generate_benchmark_plots/benchmark_105.svg
   :width: 90%
   :align: center

While not the fastest CIF parser around, ``parsnip`` achieves competetive performance when
reading files in addition to our class-leading accuracy.

Increasing Performance
^^^^^^^^^^^^^^^^^^^^^^

In some cases, particularly when constructing thousands of unit cells, the performance
of parsnip's ``build_unit_cell`` may become a bottleneck. **parsnip** includes several
tools for resolving this: first, ``parse_mode="python_float"`` attempts to build unit
cells using floating point arithmetic rather than rational expression. This is less
accurate, but is still sufficient for high-quality databases and stuctures. For the
best combination of performance and accuracy, installing the `cfractions`_ library lets
``parsnip`` use more optimized code for unit cell reconstruction. This is functionally
equivalent to the default mode, but several times faster.

.. _cfractions: https://pypi.org/project/cfractions/

.. testsetup::

    >>> import os
    >>> import numpy as np
    >>> if "doc/source" not in os.getcwd(): os.chdir("doc/source")
    >>> from parsnip import CifFile

.. doctest::

    >>> # uv pip install cfractions
    >>> from parsnip import CifFile
    >>> cif = CifFile("hP3.cif")

    >>> # If `cfractions` is available it is used by the default `parse_mode="rational"`
    >>> faster = cif.build_unit_cell(n_decimal_places=4)
    >>> faster
    array([[0.2254    , 0.        , 0.33333333],
           [0.        , 0.2254    , 0.66666667],
           [0.7746    , 0.7746    , 0.        ]])
    >>> assert faster.shape == (3, 3)



.. figure:: ../generate_benchmark_plots/benchmark_parse_modes.svg
   :width: 90%
   :align: center


Reproducing these Benchmarks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All benchmarks in this file were obtained using Python 3.13.2 on an M1 Macbook Pro, with
``parsnip`` version 0.6.1 and the ``uv.lock`` file associated with that tag. To
reproduce the results on your own hardware, run the following commands from the root
of the repository:

.. code-block:: bash

   uv sync --group tables

   # Measure parsnip's performance reading CIF files
   ./doc/generate_benchmark_plots/cif_parsing_benchmark.sh

   # Compare the efficiency of various parsing modes
   ./doc/generate_benchmark_plots/parse_mode_benchmark_plot.py

   # Measure the space group and unit cell reconstruction accuracy
   python _joss/generate_table_1.py
   python _joss/generate_table_2.py
