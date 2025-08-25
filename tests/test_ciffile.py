# ruff: noqa: SIM115
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


@pytest.mark.parametrize(
    ("input_preprocessor", "expect_warning"),
    [
        (lambda fn: open(fn), None),  # IOBase
        (lambda fn: fn, None),  # string file path
        (lambda fn: Path(fn), None),  # Path
        (lambda fn: open(fn).readlines(), None),  # list[str]
        (lambda fn: open(fn).read(), RuntimeWarning),  # raw string
    ],
)
@cif_files_mark
def test_open_methods(cif_data, input_preprocessor, expect_warning):
    print(type(input_preprocessor(cif_data.filename)))
    keys = [*cif_data.file.pairs.keys()]
    stored_data = np.asarray([*cif_data.file.pairs.values()])

    if expect_warning is not None:
        with pytest.warns(expect_warning, match="parsed as a raw CIF data block."):
            cif = CifFile(input_preprocessor(cif_data.filename))
    else:
        cif = CifFile(input_preprocessor(cif_data.filename))

    _array_assertion_verbose(keys, cif.get_from_pairs(keys), stored_data)
