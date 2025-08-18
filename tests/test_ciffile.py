import re
from pathlib import Path

import numpy as np
import pytest
from conftest import _array_assertion_verbose, cif_files_mark

from parsnip import CifFile
from parsnip._errors import ParseWarning


@cif_files_mark
def test_cast_values(cif_data):
    uncast_pairs = cif_data.file.pairs
    cif_data.file.cast_values = True

    # Casting back does nothing, but raises a warning
    expected_message = "Setting cast_values True->False has no effect on stored data."
    with pytest.warns(ParseWarning, match=expected_message):
        cif_data.file.cast_values = False

    for key, value in cif_data.file.pairs.items():
        if value == "":
            continue
        if isinstance(value, str):
            expected = uncast_pairs[key].replace("'", "").replace('"', "")
            assert re.search(r"[^0-9]|[^\.]", value) is not None
            assert value == expected
        else:
            assert isinstance(value, (int, float))

    cif_data.file._pairs = uncast_pairs  # Need to reset the data
    assert cif_data.file.pairs == uncast_pairs


@cif_files_mark
def test_open_methods(cif_data):
    keys = [*cif_data.file.pairs.keys()]
    stored_data = np.asarray([*cif_data.file.pairs.values()])

    # IOBase
    with open(cif_data.filename) as file:
        buffered = CifFile(file)
    _array_assertion_verbose(keys, buffered.get_from_pairs(keys), stored_data)

    unbuffered = CifFile(open(cif_data.filename))  # noqa: SIM115
    _array_assertion_verbose(keys, unbuffered.get_from_pairs(keys), stored_data)

    # string
    string_input = CifFile(cif_data.filename)
    _array_assertion_verbose(keys, string_input.get_from_pairs(keys), stored_data)

    # Path
    path_input = CifFile(Path(cif_data.filename))
    _array_assertion_verbose(keys, path_input.get_from_pairs(keys), stored_data)

    # Path
    path_input = CifFile(Path(cif_data.filename))
    _array_assertion_verbose(keys, path_input.get_from_pairs(keys), stored_data)

    # list[str]
    stringlist = CifFile(open(cif_data.filename).readlines())  # noqa: SIM115
    _array_assertion_verbose(keys, stringlist.get_from_pairs(keys), stored_data)
