import warnings

import numpy as np
import pytest
from ase import io
from ase.build import supercells
from conftest import box_keys, cif_files_mark
from gemmi import cif


def _gemmi_read_table(filename, keys):
    return np.array(cif.read_file(filename).sole_block().find(keys))


def _gemmi_read_keys(filename, keys, as_number=True):
    file_block = cif.read_file(filename).sole_block()
    if as_number:
        return np.array([cif.as_number(file_block.find_value(key)) for key in keys])
    return np.array([file_block.find_value(key) for key in keys])


@cif_files_mark  # TODO: test with conversions to numeric as well
def test_read_wyckoff_positions(cif_data):
    if "PDB_4INS_head.cif" in cif_data.filename:
        return
    keys = ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")
    parsnip_data = cif_data.file.read_wyckoff_positions()
    gemmi_data = _gemmi_read_table(cif_data.filename, keys)
    gemmi_data = [[cif.as_number(val) for val in row] for row in gemmi_data]
    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
def test_read_cell_params(cif_data, keys=box_keys):
    mmcif = "PDB_4INS_head.cif" in cif_data.filename
    if mmcif:
        keys = (key[0] + key[1:].replace("_", ".", 1) for key in keys)
    parsnip_data = cif_data.file.read_cell_params(mmcif=mmcif)
    gemmi_data = _gemmi_read_keys(cif_data.filename, keys)
    np.testing.assert_array_equal(parsnip_data, gemmi_data)

    if not mmcif:
        np.testing.assert_array_equal(
            cif_data.file.cell[3:], np.deg2rad(gemmi_data[3:])
        )


@cif_files_mark
def test_read_symmetry_operations(cif_data):
    if "PDB_4INS_head.cif" in cif_data.filename:
        return  # Excerpt of PDB file does not contain symmetry information

    parsnip_data = cif_data.file.read_symmetry_operations()
    gemmi_data = _gemmi_read_table(filename=cif_data.filename, keys=cif_data.symop_keys)
    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
@pytest.mark.parametrize("n_decimal_places", [3, 4, 5, 6, 9])
def test_build_unit_cell(cif_data, n_decimal_places):
    warnings.filterwarnings("ignore", "crystal system", category=UserWarning)

    if "PDB_4INS_head.cif" in cif_data.filename:
        pytest.skip("Function not compatible with PDB data.")

    parsnip_positions = cif_data.file.build_unit_cell(
        n_decimal_places=n_decimal_places, fractional=False
    )

    # Read the structure, then extract to Python builtin types. Then, wrap into the box
    ase_file = io.read(cif_data.filename)
    ase_data = supercells.make_supercell(ase_file, np.diag([1, 1, 1]))

    # Arrays must be sorted to guarantee correct comparison
    parsnip_positions = np.array(
        sorted(parsnip_positions.round(14), key=lambda x: (x[0], x[1], x[2]))
    )
    ase_positions = np.array(
        sorted(ase_data.get_positions(), key=lambda x: (x[0], x[1], x[2]))
    )

    parsnip_minmax = [parsnip_positions.min(axis=0), parsnip_positions.max(axis=0)]
    ase_minmax = [ase_positions.min(axis=0), ase_positions.max(axis=0)]
    np.testing.assert_allclose(parsnip_minmax, ase_minmax, atol=1e-12)

    np.testing.assert_allclose(parsnip_positions, ase_positions, atol=1e-12)


@cif_files_mark
def test_invalid_unit_cell(cif_data):
    if "PDB" in cif_data.filename:
        return

    previous_alpha = cif_data.file.pairs["_cell_angle_alpha"]
    cif_data.file._pairs["_cell_angle_alpha"] = "180"

    with pytest.raises(ValueError, match="outside the valid range"):
        cif_data.file.build_unit_cell()
    cif_data.file._pairs["_cell_angle_alpha"] = previous_alpha

    # cif_data.file._pairs.pop("_cell_angle_alpha")
    # with pytest.raises(ValueError, match="did not return any data"):
    #     cif_data.file.build_unit_cell()
    # cif_data.file._pairs["_cell_angle_alpha"] = previous_alpha
