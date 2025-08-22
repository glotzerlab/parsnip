import numpy as np
from conftest import (
    _gemmi_read_keys,
    all_files_mark,
    bad_cif,
    pycifrw_or_skip,
    random_keys_mark,
)
from more_itertools import flatten


def _array_assertion_verbose(keys, test_data, real_data):
    keys = np.asarray(keys)
    test_data = np.asarray(test_data)
    real_data = np.asarray(real_data)
    msg = (
        f"Key(s) {keys[test_data != real_data]} did not match:\n"
        f"{test_data[test_data != real_data]}!="
        f"{real_data[test_data != real_data]}\n"
    )
    np.testing.assert_equal(test_data, real_data, err_msg=msg)


@all_files_mark
def test_read_key_value_pairs(cif_data):
    pycif = pycifrw_or_skip(cif_data)

    invalid = [*flatten(pycif.loops.values()), *cif_data.failing]
    all_keys = [key for key in pycif.true_case.values() if key.lower() not in invalid]

    parsnip_data = cif_data.file[all_keys]
    for i, value in enumerate(parsnip_data):
        np.testing.assert_equal(cif_data.file[all_keys[i]], value)
        np.testing.assert_equal(
            cif_data.file[all_keys[i]], cif_data.file.get_from_pairs(all_keys[i])
        )
    gemmi_data = _gemmi_read_keys(cif_data.filename, keys=all_keys, as_number=False)
    _array_assertion_verbose(all_keys, parsnip_data, gemmi_data)


@all_files_mark
@random_keys_mark(n_samples=20)
def test_read_key_value_pairs_random(cif_data, keys):
    parsnip_data = np.asarray(cif_data.file[keys])
    _array_assertion_verbose(keys, parsnip_data, cif_data.file.get_from_pairs(keys))
    gemmi_data = _gemmi_read_keys(cif_data.filename, keys=keys, as_number=False)
    _array_assertion_verbose(keys, parsnip_data, gemmi_data)


def test_read_key_value_pairs_badcif(cif_data=bad_cif):
    parsnip_data = cif_data.file[cif_data.manual_keys]
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
    _array_assertion_verbose(cif_data.manual_keys, parsnip_data, correct_data)
