# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

"""Configure doctest namespace."""


import numpy as np
import pytest

from . import CifFile


# Set up doctests
@pytest.fixture(autouse=True)
def _setup_doctest(doctest_namespace):
    doctest_namespace["np"] = np
    doctest_namespace["cif"] = CifFile("doc/source/example_file.cif")
