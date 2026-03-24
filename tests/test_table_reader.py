import warnings

import numpy as np
import pytest
from ase.io import cif as asecif
from conftest import (
    _arrstrip,
    _gemmi_read_table,
    all_files_mark,
    bad_cif,
    cif_files_mark,
    pycifrw_or_skip,
)
from more_itertools import flatten

STR_WIDTH_MAX = 128
"""Maximum width for valid fields in the test suite.
Used to simplify processing of structured arrays.
"""


@all_files_mark
def test_reads_all_keys(cif_data):
    pycif = pycifrw_or_skip(cif_data)
    loop_keys = [*flatten(pycif.loops.values())]
    all_keys = [key for key in pycif.true_case.values() if key.lower() in loop_keys]

    found_labels = [*flatten(cif_data.file.loop_labels)]
    for key in all_keys:
        assert key in found_labels, f"Missing label: {key}"

    if "A2BC_tP16_76" in cif_data.filename:  # TODO: this can be supported
        pytest.xfail("Double single quote at EOL is not supported.")

    for loop in pycif.loops.values():
        loop = [pycif.true_case[key] for key in loop]
        parsnip_data = cif_data.file.get_from_loops(loop)
        gemmi_data = _gemmi_read_table(cif_data.filename, loop)
        np.testing.assert_array_equal(parsnip_data, _arrstrip(gemmi_data, r"\r"))


@cif_files_mark
def test_read_symop(cif_data):
    parsnip_data = cif_data.file.get_from_loops(cif_data.symop_keys)
    gemmi_data = _gemmi_read_table(cif_data.filename, cif_data.symop_keys)

    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@cif_files_mark
def test_read_atom_sites(cif_data):
    parsnip_data = cif_data.file.get_from_loops(cif_data.atom_site_keys)
    gemmi_data = _gemmi_read_table(cif_data.filename, cif_data.atom_site_keys)
    np.testing.assert_array_equal(parsnip_data, gemmi_data)
    assert (key in cif_data.file.loop_labels for key in cif_data.atom_site_keys)

    if "CCDC" in cif_data.filename:
        pytest.xfail("Occupancy data would need to be tiled to match ASE.")

    if cif_data.file["_atom_site_occupancy"] is not None:
        warnings.filterwarnings(
            "ignore", category=UserWarning, message="crystal system"
        )

        atoms = asecif.read_cif(cif_data.filename)

        ase_data = [
            occ for site in atoms.info["occupancy"].values() for occ in site.values()
        ]
        np.testing.assert_array_equal(
            cif_data.file.get_from_loops("_atom_site_occupancy")
            .squeeze()
            .astype(float),
            ase_data,
        )


@cif_files_mark
@pytest.mark.parametrize(
    "keys",
    [
        "*",
        "**",
        "*?*",
        "?" * 18,  # Length of '_atom_site_fract_?'
        "_space*",
        "*space*",
        "_atom_site?",
        "_atom?site*",
        "_atom_site*_?",
    ],
)
def test_wildcard_keys_loops(cif_data, keys):
    parsnip_data = cif_data.file[keys]
    raw_keys = cif_data.file._wildcard_mapping.get(keys, [keys])

    if not isinstance(parsnip_data, list):
        gemmi_data = _gemmi_read_table(cif_data.filename, raw_keys)
        np.testing.assert_array_equal(parsnip_data, gemmi_data)
        return

    start_idx = 0
    for sublist in parsnip_data:
        sublist_length = sublist.shape[1] if isinstance(sublist, np.ndarray) else 1
        keys_slice = raw_keys[start_idx : start_idx + sublist_length]
        start_idx += len(keys_slice)
        gemmi_data = _gemmi_read_table(cif_data.filename, keys_slice)
        np.testing.assert_array_equal(sublist, gemmi_data)


@cif_files_mark
@pytest.mark.parametrize(
    "subset", [[0], [1, 2, 3], [-1, 0]], ids=["single_el", "slice", "end_and_beginning"]
)
def test_partial_table_read(cif_data, subset):
    subset_of_keys = tuple(np.array(cif_data.atom_site_keys)[subset])
    parsnip_data = cif_data.file.get_from_loops(subset_of_keys)
    gemmi_data = _gemmi_read_table(cif_data.filename, subset_of_keys)

    np.testing.assert_array_equal(parsnip_data, gemmi_data)


@pytest.mark.skip("Would be nice to pass, but we are at least as good as gemmi here.")
def test_bad_cif_symop():
    parsnip_data = bad_cif.file.get_from_loops(bad_cif.symop_keys)
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
