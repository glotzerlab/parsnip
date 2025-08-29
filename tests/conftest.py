from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass
from glob import glob

import numpy as np
import pytest
from CifFile import CifFile as pycifRW
from CifFile import CifSyntaxError, StarError
from gemmi import cif

from parsnip import CifFile
from parsnip._errors import ParseWarning

ADDITIONAL_TEST_FILES_PATH = ""

rng = np.random.default_rng(seed=161181914916)

data_file_path = os.path.join(os.path.dirname(__file__), "sample_data")


def pycifrw_or_skip(cif_data):
    try:
        return pycifRW(cif_data.filename, strict=0, permissive=True).first_block()
    except StarError:
        pytest.skip("pycifRW raised a StarError!")
    except CifSyntaxError:
        pytest.skip("pycifRW raised a CifSyntaxError!")


def remove_invalid(s):
    """Our parser strips newlines and carriage returns.
    TODO: newlines should be retained
    """
    if s is None or s == "":
        return None
    return s.replace("\r", "")


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


def _value_or_nan(val):
    if val is None:
        return "_"
    return val


def _gemmi_read_keys(filename, keys, as_number=True):
    try:
        file_block = cif.read_file(filename).sole_block()
    except ValueError as e:
        if "parse error" in str(e) or "unterminated 'string'" in str(e):
            pytest.skip(f"Gemmi failed to read file: {e}")
        raise ValueError(f"Unexpected error found: {e}") from e
    except RuntimeError as e:
        if "duplicate tag" in str(e):
            pytest.skip(f"Gemmi failed to read file: {e}")
        raise RuntimeError(f"Unexpected error found: {e}") from e
    if as_number:
        return np.array(
            [cif.as_number(_value_or_nan(file_block.find_value(key))) for key in keys]
        )
    return np.array([remove_invalid(file_block.find_value(key)) for key in keys])


def _gemmi_read_table(filename, keys):
    try:
        return np.array(
            [
                [remove_invalid(x) for x in row]
                for row in cif.read_file(filename).sole_block().find(keys)
            ]
        )
    except ValueError as e:
        if "unterminated 'string'" in str(e):
            pytest.skip(f"Gemmi failed to read file: {e}")
        raise ValueError(f"Unexpected error found: {e}") from e
    except RuntimeError as e:
        if "duplicate tag" in str(e):
            pytest.skip(f"Gemmi failed to read file: {e}")
        raise RuntimeError(f"Unexpected error found: {e}") from e


def _arrstrip(arr: np.ndarray, pattern: str):
    return np.vectorize(lambda x: re.sub(pattern, "", x))(arr)


@dataclass
class CifData:
    filename: str
    file: CifFile
    symop_keys: tuple[str, ...] = ()
    atom_site_keys: tuple[str, ...] = ()
    failing: tuple[str, ...] = ()
    """Test cases that DO NOT read properly."""
    manual_keys: tuple[str, ...] = ()

    @classmethod
    def from_file(cls, file: str) -> CifData:
        """Create a CifData object from a filename and the default keys."""
        cif = cls(
            filename=os.path.join(data_file_path, file),
            file=CifFile(os.path.join(data_file_path, file)),
        )
        cif.symop_keys = (
            ("_space_group_symop_operation_xyz",)
            if cif.file["_space_group_symop_operation_xyz"] is not None
            else ("_symmetry_equiv_pos_as_xyz",)
        )
        cif.atom_site_keys = (
            *(key for key in atom_site_keys if cif.file[key] is not None),
        )
        return cif


# Assorted keys to select from
assorted_keys = np.loadtxt(os.path.join(data_file_path, "cif_file_keys.txt"), dtype=str)


def combine_marks(*marks, argnames="cif_data"):
    combinedargvalues = []
    combinedids = []
    for mark in marks:
        argvalues, ids = mark.kwargs["argvalues"], mark.kwargs["ids"]
        combinedargvalues.extend(argvalues)
        combinedids.extend(ids)
    return pytest.mark.parametrize(
        argnames=argnames,
        argvalues=combinedargvalues,
        ids=combinedids,
    )


def generate_random_key_sequences(arr, n_samples, seed=42, wildcard_probability=0):
    rng = np.random.default_rng(seed)
    wildcards = ["?", "*"]
    if wildcard_probability > 0:
        result = []
        for size in rng.integers(1, len(arr), n_samples):
            sample = rng.choice(arr, size=size, replace=False)
            for i, s in enumerate(sample):
                if rng.uniform() < wildcard_probability:
                    idx = rng.integers(0, len(s))
                    sample[i] = s[:idx] + rng.choice(wildcards) + s[idx + 1 :]
            result.append(sample)
        return result
    return [
        rng.choice(arr, size=size, replace=False)
        for size in rng.integers(1, len(arr), n_samples)
    ]


def random_keys_mark(n_samples=10):
    return pytest.mark.parametrize(
        argnames="keys",
        argvalues=generate_random_key_sequences(assorted_keys, n_samples=n_samples),
    )


def random_wildcard_keys_mark(n_samples=10):
    return pytest.mark.parametrize(
        argnames="keys",
        argvalues=generate_random_key_sequences(
            assorted_keys, n_samples=n_samples, wildcard_probability=0.9
        ),
    )


# Used for test_read_cell_params
box_keys = (
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma",
)

atom_site_keys = (
    "_atom_site_label",
    "_atom_site_type_symbol",
    "_atom_site_fract_x",
    "_atom_site_fract_y",
    "_atom_site_fract_z",
    "_atom_site_occupancy",
    "_atom_site_U_iso_or_equiv",
)

aflow_mC24 = CifData.from_file("AFLOW_mC24.cif")
amcsd_seifertite = CifData.from_file("AMCSD_meteorite.cif")
bisd_Ccmm = CifData.from_file("B-IncStrDb_Ccmm.cif")
ccdc_Pm3m = CifData.from_file("CCDC_1446529_Pm-3m.cif")
cod_aP16 = CifData.from_file("COD_1540955_aP16.cif")
cod_hP3 = CifData.from_file("COD_7228524.cif")
izasc_gismondine = CifData.from_file("zeolite_clo.cif")

pdb_4INS = CifData(
    filename=os.path.join(data_file_path, "PDB_4INS_head.cif"),
    symop_keys=("_pdbx_struct_oper_list.symmetry_operation",),
    atom_site_keys=(
        "_chem_comp.id",
        "_chem_comp.type",
        "_chem_comp.mon_nstd_flag",
        "_chem_comp.name",
        "_chem_comp.pdbx_synonyms",
        "_chem_comp.formula",
        "_chem_comp.formula_weight",
    ),
    file=CifFile(os.path.join(data_file_path, "PDB_4INS_head.cif")),
)

structure_issue_42 = CifData.from_file("no42.cif")

with pytest.warns():
    bad_cif = CifData(
        filename=os.path.join(data_file_path, "INTENTIONALLY_BAD_CIF.cif"),
        symop_keys=("_space_group_symop_id", "_space_group_symop_operation_xyz"),
        atom_site_keys=(
            "_atom_site",
            "_atom_site_type_symbol",
            "_atom_site_symmetry_multiplicity",
            "_atom_si te",
            "_atom_site_fract_z",
            "_this_key_does_not_exist",
        ),
        file=CifFile(os.path.join(data_file_path, "INTENTIONALLY_BAD_CIF.cif")),
        manual_keys=(
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
            "_cod_style_key",
            "_cif_2.0_string",
            "_cif_2.0_double_quoted_string",
            "__________asdf",
            "_-wasd",
            "not_a_valid_key",
        ),
    )
warnings.filterwarnings(
    "ignore",
    message="Duplicate key ",
    category=ParseWarning,
)
mbuild_test_files = [
    CifData.from_file(os.path.join(*fn.split(os.sep)[-2:]))
    for fn in glob(os.path.join(data_file_path, "mbuild_cifs", "*.cif"))
]


cif_data_array = [
    aflow_mC24,
    amcsd_seifertite,
    bisd_Ccmm,
    ccdc_Pm3m,
    cod_aP16,
    cod_hP3,
    izasc_gismondine,
    pdb_4INS,
    structure_issue_42,
    *mbuild_test_files,
]
cif_files_mark = pytest.mark.parametrize(
    argnames="cif_data",
    argvalues=cif_data_array,
    ids=[cif.filename.split("/")[-1] for cif in cif_data_array],
)
additional_data_array = [
    CifData(
        filename=fn,
        file=CifFile(fn),
        symop_keys=("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"),
    )
    for fn in [
        *glob(ADDITIONAL_TEST_FILES_PATH),
    ]
]
additional_files_mark = pytest.mark.parametrize(
    argnames="cif_data",
    argvalues=additional_data_array,
    ids=[cif.filename.split("/")[-1] for cif in additional_data_array],
)

all_files_mark = combine_marks(cif_files_mark, additional_files_mark)
