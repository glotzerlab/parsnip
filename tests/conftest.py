import dataclasses
import os

import numpy as np
import pytest

# ruff: noqa: N816. Allow mixed-case global variables

data_file_path = os.path.dirname(__file__) + "/sample_data/"


@dataclasses.dataclass
class CifData:
    """Class to hold the filename and stored keys for a CIF file."""

    filename: str
    symop_keys: tuple[str]
    atom_site_keys: tuple[str]
    single_value_keys: tuple[str]


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
    symop_keys=("_space_group_symop_id", "_space_group_symop_operation_xyz"),
    atom_site_keys=atom_site_keys,
    single_value_keys=(
        "_audit_creation_method",
        "_chemical_name_mineral",
        "_chemical_formula_sum",
        "_symmetry_space_group_name_H-M",
        "_aflow_title",
        "_aflow_params",
        "_aflow_params_values",
        "_aflow_Strukturbericht",
        "_aflow_Pearson",
    ),
)

bisd_Ccmm = CifData(
    filename=data_file_path + "B-IncStrDb_Ccmm.cif",
    symop_keys=("_space_group_symop_operation_xyz", "_space_group_symop_id"),
    # Our code works with extra keys, but gemmi does not!
    atom_site_keys=(atom_site_keys[0], *atom_site_keys[2:]),
    single_value_keys=(
        "_journal_name_full",
        "_journal_volume",
        "_journal_year",
        "_journal_page_first",
        "_journal_page_last",
        "_journal_paper_doi",
        "_publ_contact_author_name",
        "_publ_contact_author_email",
        "_chemical_formula_sum",
        "_space_group_crystal_system",
        "_refine_ls_wR_factor_gt",
    ),
)

ccdc_Pm3m = CifData(
    filename=data_file_path + "CCDC_1446529_Pm-3m.cif",
    symop_keys=("_space_group_symop_operation_xyz",),
    atom_site_keys=sorted(atom_site_keys),
    single_value_keys=(
        "_audit_block_doi",
        "_database_code_depnum_ccdc_archive",
        "_computing_publication_material",
        "_chemical_formula_sum",
        "_cell_formula_units_Z",
        "_space_group_crystal_system",
        "_space_group_name_H-M_alt",
        "_diffrn_ambient_temperature",
        "_reflns_number_gt",
        "_refine_ls_R_factor_gt",
        "_refine_ls_wR_factor_gt",
        "_refine_diff_density_max",
        "_refine_diff_density_min",
        "_refine_diff_density_rms",
    ),
)

cod_aP16 = CifData(
    filename=data_file_path + "COD_1540955_aP16.cif",
    symop_keys=("_symmetry_equiv_pos_as_xyz",),
    atom_site_keys=atom_site_keys,
    single_value_keys=(
        "_journal_page_first",
        "_journal_page_last",
        "_journal_volume",
        "_journal_year",
        "_chemical_formula_sum",
        "_chemical_name_systematic",
        "_space_group_IT_number",
        "_symmetry_space_group_name_Hall",
        "_symmetry_space_group_name_H-M",
        "_cell_formula_units_Z",
        "_cell_volume",
        "_citation_journal_id_ASTM",
        "_cod_data_source_file",
        "_cod_data_source_block",
        "_cod_original_cell_volume",
        "_cod_original_formula_sum",
        "_cod_database_code",
    ),
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
    single_value_keys=(
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

cif_data_array = [aflow_mC24, bisd_Ccmm, ccdc_Pm3m, cod_aP16]
cif_files_mark = pytest.mark.parametrize(
    argnames="cif_data",
    argvalues=cif_data_array,
    ids=[cif.filename.split("/")[-1] for cif in cif_data_array],
)
