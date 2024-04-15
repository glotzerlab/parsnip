import numpy as np
from conftest import box_keys, cif_files_mark
from gemmi import cif

from parsnip.unitcells import (
    read_cell_params,
    read_fractional_positions,
)


def _gemmi_read_table(filename, keys):
    return np.array(cif.read_file(filename).sole_block().find(keys))


def _gemmi_read_keys(filename, keys, as_number=True):
    file_block = cif.read_file(filename).sole_block()
    if as_number:
        return np.array([cif.as_number(file_block.find_value(key)) for key in keys])
    else:
        return np.array([file_block.find_value(key) for key in keys])


@cif_files_mark
def test_read_fractional_positions(cif_data):
    if "PDB_4INS_head.cif" in cif_data.filename:
        return
    keys = ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")
    parsnip_data = read_fractional_positions(filename=cif_data.filename)
    gemmi_data = _gemmi_read_table(cif_data.filename, keys)
    gemmi_data = [[cif.as_number(val) for val in row] for row in gemmi_data]
    np.testing.assert_allclose(parsnip_data, gemmi_data)


@cif_files_mark
def test_read_cell_params(cif_data, keys=box_keys):
    mmcif = "PDB_4INS_head.cif" in cif_data.filename
    parsnip_data = read_cell_params(filename=cif_data.filename, mmcif=mmcif)
    if mmcif:
        keys = (key[0] + key[1:].replace("_", ".", 1) for key in keys)
    gemmi_data = _gemmi_read_keys(cif_data.filename, keys)
    np.testing.assert_array_equal(parsnip_data, gemmi_data)
