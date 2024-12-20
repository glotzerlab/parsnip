"""``parsnip``: a package for the simple reading and processing of .cif files.

While there are many packages for handling cif files exist, the vast majority suffer
from decades of feature creep and high levels of complexity. ``parsnip`` provides a
simple and minimal interface for reading cif files into Python primitive data structures
and numpy arrays. The ``parsnip.parse`` module contains exactly two functions that read
key-value and tabular data from cif files, and are all that are required for most users.
The ``parsnip.patterns`` module includes a few convience features for manipulation of
the read data, and the ``parsnip.unitcells`` module includes functions to reconstruct a
crystal's unit cell's basis positions from data stored in cif files.
"""
from . import parse, patterns, unitcells

__version__ = "0.0.2"
