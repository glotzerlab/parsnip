import numpy as np
import pytest
from conftest import bad_cif, cif_files_mark, random_keys_mark
from gemmi import cif
from ase.io import cif as asecif


def _gemmi_read_keys(filename, keys, as_number=True):
    file_block = cif.read_file(filename).sole_block()
    if as_number:
        return np.array([cif.as_number(file_block.find_value(key)) for key in keys])
    else:
        return np.array([file_block.find_value(key) for key in keys])


@cif_files_mark
def test_read_key_value_pairs(cif_data):
    parsnip_data = cif_data.file[cif_data.single_value_keys]
    gemmi_data = _gemmi_read_keys(
        cif_data.filename, keys=cif_data.single_value_keys, as_number=False
    )
    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
@random_keys_mark(n_samples=20)
def test_read_key_value_pairs_random(cif_data, keys):
    parsnip_data = cif_data.file[keys]
    gemmi_data = _gemmi_read_keys(cif_data.filename, keys=keys, as_number=False)
    np.testing.assert_array_equal(parsnip_data, gemmi_data)


def test_read_key_value_pairs_badcif(cif_data=bad_cif):
    parsnip_data = cif_data.file[cif_data.single_value_keys]
    correct_data = [
        "1.000000(x)",
        "4.32343242",
        "3.1415926535897932384626433832795028841971693993751058209749",
        "90.00000",
        "-10.12345",
        "210.00000",
        "\\t _1.234-56789",
        r"45.6a/\s",
        None,
    ]
    np.testing.assert_array_equal(parsnip_data, correct_data)
