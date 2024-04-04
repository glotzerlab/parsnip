import dataclasses
import os

import pytest

# ruff: noqa: N816. Allow mixed-case global variables


@dataclasses.dataclass
class CifData:
    """Class to hold the filename and stored keys for a CIF file."""

    filename: str
    symop_keys: tuple[str]
    atom_site_keys: tuple[str]


box_keys = (
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma",
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
)
atom_site_keys = (
    "_atom_site_label",
    "_atom_site_type_symbol",
    "_atom_site_fract_x",
    "_atom_site_fract_y",
    "_atom_site_fract_z",
    "_atom_site_occupancy",
)


data_file_path = os.path.dirname(__file__) + "/sample_data/"


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
    ),
)

cif_data_array = [aflow_mC24, bisd_Ccmm, ccdc_Pm3m, cod_aP16]
cif_files_mark = pytest.mark.parametrize(
    argnames="cif_data",
    argvalues=cif_data_array,
    ids=[cif.filename.split("/")[-1] for cif in cif_data_array],
)
