import re

import numpy as np

from parsnip import CifFile
from parsnip._errors import ParseWarning
from conftest import cif_files_mark

import pytest

from copy import deepcopy

# @pytest.mark.skip("cast_values changes the data for other tests.")
@cif_files_mark
def test_cast_values(cif_data):
    uncast_pairs = deepcopy(cif_data.file.pairs)
    cif_data.file.cast_values = True

    # Casting back does nothing, but raises a warning
    expected_message = "Setting cast_values True->False has no effect on stored data."
    with pytest.warns(ParseWarning, match=expected_message):
        cif_data.file.cast_values = False

    for key, value in cif_data.file.pairs.items():
        if value == "":
            continue
        if isinstance(value, str):
            expected = deepcopy(uncast_pairs[key]).replace("'", "").replace('"', "")
            assert re.search(r"[^0-9]|[^\.]", value) is not None
            assert value == expected
        else:
            assert isinstance(value, (int, float))

    cif_data.file._pairs = uncast_pairs # Need to reset the data
    assert cif_data.file.pairs == uncast_pairs
