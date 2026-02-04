import json
from pathlib import Path

import numpy as np
import pytest

import parsnip
from parsnip import CifFile
from parsnip.patterns import _normalize

# Load reference data, relative to the main package installation
SYMOPS_PATH = Path(parsnip.__file__).parent / "symops.json"
with open(SYMOPS_PATH) as f:
    RAW_DATA = json.load(f)


def lookup_test_cases():
    """Construct all possible space group queries from `symops.json`."""
    cases = []

    # Hall Symbols
    for hall, data in RAW_DATA.items():
        cases.append(("_space_group_name_Hall", hall, data["symops"]))

    # H-M Symbols, tiled out into the full/short and with/without setting options.
    hm_map = {}
    for data in RAW_DATA.values():
        variants = [
            data["hermann_mauguin_full"],
            data["hermann_mauguin_full"].split(":")[0],
            data["hermann_mauguin_short"],
            data["hermann_mauguin_short"].split(":")[0],
        ]
        for v in variants:
            hm_map[_normalize(v)] = (v, data["symops"])

    for _, (original_hm, symops) in hm_map.items():
        cases.append(("_space_group_name_H-M_alt", original_hm, symops))
        cases.append(("_symmetry_space_group_name_H-M", original_hm, symops))

    # International tables numbers
    it_map = {}
    for data in RAW_DATA.values():
        it_map[_normalize(data["table_number"])] = (
            data["table_number"],
            data["symops"],
        )

    for _, (original_it, symops) in it_map.items():
        cases.append(("_space_group_IT_number", original_it, symops))
        cases.append(("_symmetry_Int_Tables_number", original_it, symops))

    return cases


@pytest.mark.parametrize(("key", "value", "expected_symops"), lookup_test_cases())
def test_symops_lookup(key, value, expected_symops):
    # Construct a simple CIF with the key-value pair. This is parsed as a raw string.
    cif_content = f"data_test\n{key} '{value}'"

    with pytest.warns(RuntimeWarning, match="File input was parsed"):
        cif = CifFile(cif_content)

    assert cif.symops is not None
    generated_symops = cif.symops.flatten().tolist()
    np.testing.assert_array_equal(generated_symops, expected_symops)


def test_symops_lookup_none():
    cif_content = "data_test\n_cell_length_a 10\nloop_\n_atom_site_label\nC"
    with pytest.warns(RuntimeWarning, match="File input was parsed"):
        cif = CifFile(cif_content)
    assert cif.symops is None
