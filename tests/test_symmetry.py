import numpy as np
import pytest

from parsnip import CifFile


def test_symops_lookup_it_number():
    cif_content = """
data_
_space_group_IT_number 2
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C 0 0 0
    """
    with pytest.warns(RuntimeWarning, match="File input was parsed"):
        cif = CifFile(cif_content)
    symops = cif.symops
    assert symops is not None
    assert isinstance(symops, np.ndarray)
    assert symops.shape == (2, 1)
    expected = ["x,y,z", "-x,-y,-z"]
    flattened = symops.flatten()
    assert flattened[0] == expected[0]
    assert flattened[1] == expected[1]


def test_symops_lookup_int_tables_number():
    cif_content = """
data_test
_symmetry_Int_Tables_number 2
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C 0 0 0
    """
    with pytest.warns(RuntimeWarning, match="File input was parsed"):
        cif = CifFile(cif_content)
    symops = cif.symops
    assert symops is not None
    assert symops.shape == (2, 1)
    expected = ["x,y,z", "-x,-y,-z"]
    flattened = symops.flatten()
    assert flattened[0] == expected[0]
    assert flattened[1] == expected[1]


def test_symops_lookup_none():
    cif_content = """
data_test
# No group number
_cell_length_a 10
loop_
_atom_site_label
C
    """
    with pytest.warns(RuntimeWarning, match="File input was parsed"):
        cif = CifFile(cif_content)
    assert cif.symops is None
