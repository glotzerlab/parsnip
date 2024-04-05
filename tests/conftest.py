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
)

bisd_Ccmm = CifData(
    filename=data_file_path + "B-IncStrDb_Ccmm.cif",
    symop_keys=("_space_group_symop_operation_xyz", "_space_group_symop_id"),
    # Our code works with extra keys, but gemmi does not!
    atom_site_keys=(atom_site_keys[0], *atom_site_keys[2:]),
)

ccdc_Pm3m = CifData(
    filename=data_file_path + "CCDC_1446529_Pm-3m.cif",
    symop_keys=("_space_group_symop_operation_xyz",),
    atom_site_keys=sorted(atom_site_keys),
)

cod_aP16 = CifData(
    filename=data_file_path + "COD_1540955_aP16.cif",
    symop_keys=("_symmetry_equiv_pos_as_xyz",),
    atom_site_keys=atom_site_keys,
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
)

cif_data_array = [aflow_mC24, bisd_Ccmm, ccdc_Pm3m, cod_aP16]
cif_files_mark = pytest.mark.parametrize(
    argnames="cif_data",
    argvalues=cif_data_array,
    ids=[cif.filename.split("/")[-1] for cif in cif_data_array],
)
