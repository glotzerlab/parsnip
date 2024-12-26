import warnings

import numpy as np
import pytest
from ase.io import cif as asecif
from conftest import bad_cif, cif_files_mark
from gemmi import cif

STR_WIDTH_MAX = 128
"""Maximum width for valid fields in the test suite.
Used to simplify processing of structured arrays.
"""


def _gemmi_read_table(filename, keys):
    return np.array(cif.read_file(filename).sole_block().find(keys))


@cif_files_mark
def test_read_symop(cif_data):
    parsnip_data = cif_data.file.get_from_tables(cif_data.symop_keys)
    gemmi_data = _gemmi_read_table(cif_data.filename, cif_data.symop_keys)

    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
def test_read_atom_sites(cif_data):
    parsnip_data = cif_data.file.get_from_tables(cif_data.atom_site_keys)
    gemmi_data = _gemmi_read_table(cif_data.filename, cif_data.atom_site_keys)
    if "CCDC" not in cif_data.filename and "PDB" not in cif_data.filename:
        # These CCDC and PDB files cannot be read by ASE
        warnings.filterwarnings("ignore", category=UserWarning)
        atoms = asecif.read_cif(cif_data.filename)
        ase_data = [
            occ for site in atoms.info["occupancy"].values() for occ in site.values()
        ]
        np.testing.assert_array_equal(
            cif_data.file.get_from_tables("_atom_site_occupancy")
            .squeeze()
            .astype(float),
            ase_data,
        )
    np.testing.assert_array_equal(parsnip_data, gemmi_data)
    assert (key in cif_data.file.table_labels for key in cif_data.atom_site_keys)


@cif_files_mark
@pytest.mark.parametrize(
    "subset", [[0], [1, 2, 3], [4, 0]], ids=["single_el", "slice", "end_and_beginning"]
)
def test_partial_table_read(cif_data, subset):
    subset_of_keys = tuple(np.array(cif_data.atom_site_keys)[subset])
    parsnip_data = cif_data.file.get_from_tables(subset_of_keys)
    gemmi_data = _gemmi_read_table(cif_data.filename, subset_of_keys)

    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@pytest.mark.skip("Would be nice to pass, but we are at least as good as gemmi here.")
def test_bad_cif_symop(cif_data=bad_cif):
    # This file is thouroughly cooked - gemmi will not even read it.
    parsnip_data = cif_data.file.get_from_tables(cif_data.symop_keys)
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


@pytest.mark.skip("Too corrupted to be read")
def test_bad_cif_atom_sites(cif_data=bad_cif):
    parsnip_data = cif_data.file[cif_data.atom_site_keys]
    np.testing.assert_array_equal(
        parsnip_data[:, 0],
        np.array(["Aa(3)", "SL", "Oo", "O0f"]),
    )
    # "_atom_site_type_symbol"
    np.testing.assert_array_equal(parsnip_data[:, 1], ["Bb", "SM", "O", "O"])

    # "_atom_site_symmetry_multiplicity"
    np.testing.assert_array_equal(parsnip_data[:, 2], ["1", "3", "5", "7"])

    # "_atom_si te"
    np.testing.assert_array_equal(
        parsnip_data[:, 3], ["0.00000(1)", "0.00000", "0.19180", "0.09390"]
    )

    # "_atom_site_fract_z"
    np.testing.assert_array_equal(
        parsnip_data[:, 4], ["0.25000", "0.(28510)", "0.05170", "0.41220"]
    )
