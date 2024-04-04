import numpy as np
import pytest
from conftest import bad_cif, cif_files_mark
from gemmi import cif

from parsnip.parse import read_table


def _gemmi_read_table(filename, keys):
    return np.array(cif.read_file(filename).sole_block().find(keys))


@cif_files_mark
def test_read_symop(cif_data):
    parsnip_data = read_table(filename=cif_data.filename, keys=cif_data.symop_keys)
    gemmi_data = _gemmi_read_table(cif_data.filename, cif_data.symop_keys)

    # We replace ", " strings with "," to ensure data is collected properly
    # We have to apply this same transformation to the gemmi data to check correctness.
    if "CCDC_1446529_Pm-3m.cif" in cif_data.filename:
        gemmi_data = np.array(
            [[item.replace(", ", ",") for item in row] for row in gemmi_data]
        )

    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
def test_read_atom_sites(cif_data):
    parsnip_data = read_table(
        filename=cif_data.filename,
        keys=cif_data.atom_site_keys,
    )
    gemmi_data = _gemmi_read_table(cif_data.filename, cif_data.atom_site_keys)

    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
@pytest.mark.parametrize(
    "subset", [[0], [1, 2, 3], [4, 0]], ids=["single_el", "slice", "end_and_beginning"]
)
def test_partial_table_read(cif_data, subset):
    subset_of_keys = tuple(np.array(cif_data.atom_site_keys)[subset])
    parsnip_data = read_table(
        filename=cif_data.filename,
        keys=subset_of_keys,
    )
    gemmi_data = _gemmi_read_table(cif_data.filename, subset_of_keys)

    np.testing.assert_array_equal(parsnip_data, gemmi_data)


def test_bad_cif_symop(cif_data=bad_cif):
    # This file is thouroughly cooked - gemmi will not even read it.
    parsnip_data = read_table(
        filename=cif_data.filename,
        keys=cif_data.symop_keys,
    )
    correct_data = [
        ["1", "x,y,z"],
        ["2", "-x,y,-z*1/2"],
        ["3", "-x,-y,-z"],
        ["4", "x,=y,z/1/2"],
        ["5", "x-1/2,y+1/2,z"],
        ["6", "-x+1/2,ya1/2,-z+1/2"],
        ["7", "-x+1/2,-y81/2,-z"],
        ["8", "x+1/2,-y+1/2,z01/2"],
    ]

    np.testing.assert_array_equal(parsnip_data, correct_data)


def test_bad_cif_atom_sites(cif_data=bad_cif):
    parsnip_data = read_table(
        filename=cif_data.filename,
        keys=cif_data.atom_site_keys,
    )
    # This file is thouroughly cooked - gemmi will not even read it.
    print(parsnip_data)
