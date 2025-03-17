import re
import warnings
from contextlib import nullcontext
from glob import glob

import numpy as np
import pytest
from ase import io
from ase.build import supercells
from conftest import (
    ADDITIONAL_TEST_FILES_PATH,
    _arrstrip,
    _gemmi_read_keys,
    all_files_mark,
    box_keys,
    cif_files_mark,
)
from gemmi import cif
from more_itertools import flatten

from parsnip import CifFile


def _gemmi_read_table(filename, keys):
    try:
        return np.array(cif.read_file(filename).sole_block().find(keys))
    except (RuntimeError, ValueError):
        pytest.xfail("Gemmi failed to read file!")


@all_files_mark  # TODO: test with conversions to numeric as well
def test_read_wyckoff_positions(cif_data):
    if "PDB_4INS_head.cif" in cif_data.filename:
        return
    keys = ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")
    parsnip_data = cif_data.file.wyckoff_positions
    gemmi_data = _gemmi_read_table(cif_data.filename, keys)
    gemmi_data = [[cif.as_number(val) for val in row] for row in gemmi_data]
    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@all_files_mark
def test_read_cell_params(cif_data, keys=box_keys):
    if "PDB_4INS_head.cif" in cif_data.filename:
        keys = (key[0] + key[1:].replace("_", ".", 1) for key in keys)
    parsnip_data = cif_data.file.read_cell_params()
    gemmi_data = _gemmi_read_keys(cif_data.filename, keys)
    np.testing.assert_array_equal(parsnip_data, gemmi_data)

    normalized = cif_data.file.read_cell_params(normalize=True)
    assert normalized[3:] == parsnip_data[3:]  # Should not change the angles
    assert min(normalized[:3]) == 1


@all_files_mark
def test_read_symmetry_operations(cif_data):
    if "PDB_4INS_head.cif" in cif_data.filename:
        return  # Excerpt of PDB file does not contain symmetry information

    parsnip_data = cif_data.file.symops
    gemmi_data = [
        _gemmi_read_table(filename=cif_data.filename, keys=[k])
        for k in cif_data.symop_keys
    ]
    gemmi_data = gemmi_data[0] if len(gemmi_data[0]) != 0 else gemmi_data[1]
    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
@pytest.mark.parametrize("n_decimal_places", [3, 4, 6, 9])
@pytest.mark.parametrize(
    "cols",
    [
        None,
        "_atom_site_type_symbol",
        ["_atom_site_type_symbol", "_atom_site_occupancy"],
    ],
)
def test_build_unit_cell(cif_data, n_decimal_places, cols):
    warnings.filterwarnings("ignore", "crystal system", category=UserWarning)

    if "PDB_4INS_head.cif" in cif_data.filename:
        return

    should_raise = cols is not None and any(
        k not in flatten(cif_data.file.loop_labels) for k in np.atleast_1d(cols)
    )
    occupancies, read_data = None, None
    with (
        pytest.raises(ValueError, match=r"not included in the `_atom_site_fract_\[xyz")
        if should_raise
        else nullcontext()
    ):
        read_data = cif_data.file.build_unit_cell(
            n_decimal_places=n_decimal_places, additional_columns=cols
        )

    if read_data is None:
        return  # ValueError was raised

    if cols is None:
        parsnip_positions = read_data @ cif_data.file.lattice_vectors.T
    else:
        auxiliary_arr, parsnip_positions = read_data
        parsnip_positions = parsnip_positions @ cif_data.file.lattice_vectors.T

        che_symbols = _arrstrip(auxiliary_arr[:, 0], r"[^A-Za-z]+")
        if isinstance(cols, list):
            occupancies = _arrstrip(auxiliary_arr[:, 1], r"[^\d\.]+").astype(float)

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
    ase_symbols = np.array(ase_data.get_chemical_symbols())

    parsnip_minmax = [parsnip_positions.min(axis=0), parsnip_positions.max(axis=0)]
    ase_minmax = [ase_positions.min(axis=0), ase_positions.max(axis=0)]
    np.testing.assert_allclose(parsnip_minmax, ase_minmax, atol=1e-12)

    if cols is not None:
        # NOTE: ASE saves the occupancies of the most dominant species!
        # Parsnip makes no assumptions regarding the correct occupancy
        # Check all if full occupancy, partial if occ is ndarray and None if occ is None
        mask = (
            ... if cif_data.file["_atom_site_occupancy"] is None else occupancies == 1
        )
        np.testing.assert_equal(che_symbols[mask], ase_symbols[mask])

    if "zeolite" in cif_data.filename:
        return  # Four decimal places not sufficient to reconstruct this structure
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


@pytest.mark.skipif(ADDITIONAL_TEST_FILES_PATH == "", reason="No test path provided.")
@pytest.mark.parametrize("filename", glob(ADDITIONAL_TEST_FILES_PATH))
def test_build_accuracy(filename):
    def n_from_pearson(p: str) -> int:
        return int(re.sub(r"[^\w]", "", p)[2:])

    cif = CifFile(filename)
    n, uc = n_from_pearson(cif["*Pearson"]), cif.build_unit_cell()
    uc = np.array(sorted(uc, key=lambda x: tuple(x)))
    # np.testing.assert_equal(
    #     len(uc), n, err_msg="cell does not match Pearson symbol!"
    # )
