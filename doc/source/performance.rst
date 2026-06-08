.. _performance:

===========
Performance
===========

Although `parsnip` is optimized for accuracy rather than performace, our parsing
strategy results in file reads and unit-cell reconstructions that are 3-500x faster than
comparable tools. The following data makes use of [ltalirz/cif-parsing-benchmark](https://github.com/ltalirz/cif-parsing-benchmark), a benchmark that profiles the parsing
throughput of several open-source CIF libraries.

.. figure:: ../generate_benchmark_plots/benchmark_105.svg
   :width: 50%
   :align: center
