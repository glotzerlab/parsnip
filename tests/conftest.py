# ruff: noqa: N816. Allow mixed-case global variables
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pytest

from parsnip import CifFile

rng = np.random.default_rng(seed=161181914916)

data_file_path = os.path.dirname(__file__) + "/sample_data/"


@dataclass
class CifData:
    filename: str
    file: CifFile
    symop_keys: tuple[str, ...] = ()
    atom_site_keys: tuple[str, ...] = ()
    failing: tuple[str, ...] = ()
    """Test cases that DO NOT read properly."""
    manual_keys: tuple[str, ...] = ()


# Assorted keys to select from
assorted_keys = np.loadtxt(data_file_path + "cif_file_keys.txt", dtype=str)


def generate_random_key_sequences(arr, n_samples, seed=42):
    rng = np.random.default_rng(seed)
    return [
        rng.choice(arr, size=size, replace=False)
        for size in rng.integers(1, len(arr), n_samples)
    ]


def random_keys_mark(n_samples=10):
    return pytest.mark.parametrize(
        argnames="keys",
        argvalues=generate_random_key_sequences(assorted_keys, n_samples=n_samples),
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
)

aflow_mC24 = CifData(
    filename=data_file_path + "AFLOW_mC24.cif",
    symop_keys=("_space_group_symop_operation_xyz",),
    atom_site_keys=atom_site_keys,
    file=CifFile(data_file_path + "AFLOW_mC24.cif"),
)

amcsd_seifertite = CifData(
    filename=data_file_path + "AMCSD_meteorite.cif",
    symop_keys=("_space_group_symop_operation_xyz",),
    atom_site_keys=(
        atom_site_keys[0],
        *atom_site_keys[2:5],
        "_atom_site_U_iso_or_equiv",
    ),
    file=CifFile(data_file_path + "AMCSD_meteorite.cif"),
)

bisd_Ccmm = CifData(
    filename=data_file_path + "B-IncStrDb_Ccmm.cif",
    symop_keys=("_space_group_symop_operation_xyz",),
    atom_site_keys=(atom_site_keys[0], atom_site_keys[-1], *atom_site_keys[2:-1]),
    file=CifFile(data_file_path + "B-IncStrDb_Ccmm.cif"),
)

ccdc_Pm3m = CifData(
    filename=data_file_path + "CCDC_1446529_Pm-3m.cif",
    symop_keys=("_space_group_symop_operation_xyz",),
    atom_site_keys=(*sorted(atom_site_keys),),
    file=CifFile(data_file_path + "CCDC_1446529_Pm-3m.cif"),
    failing=("_refine_ls_weighting_details",),
)

cod_aP16 = CifData(
    filename=data_file_path + "COD_1540955_aP16.cif",
    symop_keys=("_symmetry_equiv_pos_as_xyz",),
    atom_site_keys=atom_site_keys,
    file=CifFile(data_file_path + "COD_1540955_aP16.cif"),
    failing=("_journal_name_full",),
)

izasc_gismondine = CifData(
    filename=data_file_path + "zeolite_clo.cif",
    symop_keys=("_symmetry_equiv_pos_as_xyz",),
    atom_site_keys=atom_site_keys[:-1],
    file=CifFile(data_file_path + "zeolite_clo.cif"),
)

# with pytest.warns(ParseWarning, match="cannot be resolved into a table"):
pdb_4INS = CifData(
    filename=data_file_path + "PDB_4INS_head.cif",
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
    file=CifFile(data_file_path + "PDB_4INS_head.cif"),
)

bad_cif = CifData(
    filename=data_file_path + "INTENTIONALLY_BAD_CIF.cif",
    symop_keys=("_space_group_symop_id", "_space_group_symop_operation_xyz"),
    atom_site_keys=(
        "_atom_site",
        "_atom_site_type_symbol",
        "_atom_site_symmetry_multiplicity",
        "_atom_si te",
        "_atom_site_fract_z",
        "_this_key_does_not_exist",
    ),
    file=CifFile(data_file_path + "INTENTIONALLY_BAD_CIF.cif"),
    manual_keys=(
        "_cell_length_a",
        "_cell_length_b",
        "_cell_length_c",
        "_cell_angle_alpha",
        "_cell_angle_beta",
        "_cell_angle_gamma",
        "__________asdf",
        "_-wasd",
        "not_a_valid_key",
    ),
)

cif_data_array = [
    aflow_mC24,
    amcsd_seifertite,
    bisd_Ccmm,
    ccdc_Pm3m,
    cod_aP16,
    izasc_gismondine,
    pdb_4INS,
]
cif_files_mark = pytest.mark.parametrize(
    argnames="cif_data",
    argvalues=cif_data_array,
    ids=[cif.filename.split("/")[-1] for cif in cif_data_array],
)
LINE_TEST_CASES = [
    None,
    "_key",
    "__key",
    "_key.loop_",
    "asdf        ",
    "loop_",
    "",
    " ",
    "# comment",
    "_key#comment_ # 2",
    "loop_##com",
    "'my quote' # abc",
    "\"malformed ''#",
    ";oddness\"'\n;asdf",
    "_key.loop.inner_",
    "normal_case",
    "multi.periods....",
    "__underscore__",
    "_key_with_numbers123",
    "test#hash",
    "#standalone",
    "'quote_in_single'",
    '"quote_in_double"',
    " \"mismatched_quotes' ",
    ";semicolon_in_text",
    ";;double_semicolon",
    "trailing_space ",
    " leading_space",
    "_key.with#hash.loop",
    "__double#hash#inside__",
    "single'; quote",
    'double;"quote',
    "#comment;inside",
    "_tricky'combination;#",
    ";'#another#combo;'",
    '"#edge_case"',
    'loop;"#complex"case',
    "_'_weird_key'_",
    "semi;;colon_and_hash#",
    "_odd.key_with#hash",
    "__leading_double_underscore__",
    "middle;;semicolon",
    "#just_a_comment",
    '"escaped "quote"',
    "'single_quote_with_hash#'",
    "_period_end.",
    "loop_.trailing_",
    "escaped\\nnewline",
    "#escaped\\twith_tab ",
    "only;semicolon",
    "trailing_semicolon;",
    "leading_semicolon;",
    "_key;.semicolon",
    "semicolon;hash#",
    "complex\"';hash#loop",
    "just_text",
    'loop#weird"text;',
    "nested'quotes\"here   ",
    "normal_case2",
    "__underscored_case__",
    'escaped\\"quotes#',
    ";semicolon#hash;",
    "double#hash_inside##key",
    "__double..periods__",
    "key#comment ; and_more",
    "_weird_;;#combo    ",
]
