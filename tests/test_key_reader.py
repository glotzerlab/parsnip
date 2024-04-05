import numpy as np
from gemmi import cif

from conftest import box_keys, cif_files_mark, random_keys_mark
from parsnip.parse import read_cell_params, read_key_value_pairs


def _gemmi_read_keys(filename, keys, as_number=True):
    file_block = cif.read_file(filename).sole_block()
    if as_number:
        return np.array([cif.as_number(file_block.find_value(key)) for key in keys])
    else:
        return np.array([file_block.find_value(key) for key in keys])


@cif_files_mark
def test_read_cell_params(cif_data, keys=box_keys):
    parsnip_data = read_cell_params(filename=cif_data.filename)
    gemmi_data = _gemmi_read_keys(cif_data.filename, keys)
    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
@random_keys_mark(n_samples=20)
def test_read_cell_params(cif_data, keys):
    parsnip_data = read_key_value_pairs(filename=cif_data.filename, keys=keys)
    gemmi_data = _gemmi_read_keys(cif_data.filename, keys=keys, as_number=False)
    np.testing.assert_array_equal([*parsnip_data.values()], gemmi_data)
