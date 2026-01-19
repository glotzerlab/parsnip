# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

"""Configure doctest namespace."""

import numpy as np
import pytest

from . import CifFile


# Set up doctests
@pytest.fixture(autouse=True)
def _setup_doctest(doctest_namespace):
    import os

    if "doc/source" not in os.getcwd():
        os.chdir("doc/source")
    doctest_namespace["np"] = np
    doctest_namespace["cif"] = CifFile("example_file.cif")
