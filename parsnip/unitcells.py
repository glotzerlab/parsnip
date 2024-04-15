"""Functions for constructing unit cells from cif file data.

Rather than storing an entire unit cell's atomic positions, cif files instead include
the data required to recreate those positions based on symmetry rules. Symmetry
operations (stored as strings of x,y,z position permutations) are applied to the Wyckoff
(symmetry irreducible) positions to create a list of possible atomic sites. These are
then wrapped into the unit cell and filtered for uniqueness to yield the final crystal.

"""
from __future__ import annotations

import re
import warnings

import numpy as np

from parsnip._errors import ParseWarning
from parsnip._utils import _deg2rad
from parsnip.parse import read_key_value_pairs, read_table
from parsnip.patterns import cast_array_to_float


def read_fractional_positions(
    filename: str,
    regex_filter: tuple[tuple[str, str]] | None = ((r",\s+", ",")),
):
    r"""Extract the fractional X,Y,Z coordinates from a CIF file.

    .. warning::

        This function ONLY returns the symmetry irreducible positions that are directly
        stored in the CIF file. To build out the full unit cell, use
        :meth:`extract_unit_cell`.

    Args:
        filename (str): The name of the .cif file to be parsed.
        regex_filter (tuple[tuple[str]], optional):
            A tuple of strings that are compiled to a regex filter and applied to each
            data line. Default value = ``None``

    Returns:
        :math:`(N, 3)` :class:`numpy.ndarray[np.float32]`:
            Fractional X,Y,Z coordinates of the unit cell.
    """
    xyz_keys = ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")
    xyz_data = read_table(
        filename=filename,
        keys=xyz_keys,
        nondelimiting_whitespace_replacement="",
        regex_filter=regex_filter,
    )
    xyz_data = cast_array_to_float(arr=xyz_data, dtype=np.float64)
    assert xyz_data.shape[1] == 3

    return xyz_data


def read_symmetry_operations(
    filename,
    regex_filter: tuple[tuple[str, str]] | None = None,
):
    r"""Extract the symmetry operations from a CIF file.

    Args:
        filename (str): The name of the .cif file to be parsed.
        regex_filter (tuple[tuple[str]], optional):
            A tuple of strings that are compiled to a regex filter and applied to each
            data line. Default value = ``None``

    Returns:
        :math:`(N,)` :class:`numpy.ndarray[str]`:
            Symmetry operations as strings.
    """
    symmetry_keys = (
        "_symmetry_equiv_pos_as_xyz",
        "_space_group_symop_operation_xyz",
    )
    with warnings.catch_warnings(category=ParseWarning, action="ignore"):
        # Only one of the two keys will be matched. We can safely ignore that warning.
        data = read_table(
            filename=filename,
            keys=symmetry_keys,
            regex_filter=regex_filter,
            nondelimiting_whitespace_replacement="",
        )

    return data


def read_cell_params(filename, degrees: bool = True, mmcif: bool = False):
    r"""Read the cell lengths and angles from a CIF file.

    Args:
        filename (str): The name of the .cif file to be parsed.
        degrees (bool, optional):
            When True, angles are returned in degrees (as per the cif spec). When False,
            angles are converted to radians.
            Default value = ``True``
        mmcif (bool, optional):
            When False, the standard CIF key naming is used (e.g. _cell_angle_alpha).
            When True, the mmCIF standard is used instead (e.g. cell.angle_alpha).
            Default value = ``False``

    Returns:
        tuple:
            The box vector lengths and angles in degrees or radians
            :math:`(L_1, L_2, L_3, \alpha, \beta, \gamma)`.
    """
    if mmcif:
        angle_keys = ("_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma")
        box_keys = ("_cell.length_a", "_cell.length_b", "_cell.length_c") + angle_keys
    else:
        angle_keys = ("_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")
        box_keys = ("_cell_length_a", "_cell_length_b", "_cell_length_c") + angle_keys
    cell_data = read_key_value_pairs(filename, keys=box_keys, only_read_numerics=True)

    assert all(value is not None for value in cell_data.values())
    assert all(0 < cell_data[key] < 180 for key in angle_keys)

    if not degrees:
        for key in angle_keys:
            cell_data[key] = _deg2rad(cell_data[key])

    return tuple(cell_data.values())


def _safe_eval(str_input: str, x: int | float, y: int | float, z: int | float):
    """Attempt to safely evaluate a string of symmetry equivalent positions.

    Python's ``eval`` is notoriously unsafe. While we could evaluate the entire list at
    once, doing so carries some risk. The typical alternative, ``ast.literal_eval``,
    doesnot work because we need to evaluate mathematical operations.

    We first replace the x,y,z values with ordered fstring inputs, to simplify the input
    of fractional coordinate data. This is done for convenience more than security.

    Once we substitute in the x,y,z values, we should have a string version of a list
    containing only numerics and math operators. We apply a substitution to ensure this
    is the case, then perform one final check. If it passes, we evaluate the list. Note
    that __builtins__ is set to {}, meaning importing functions is not possible. The
    __locals__ dict is also set to {}, so no variables are accessible in the evaluation.

    I cannot guarantee this is fully safe, but it at the very least makes it extremely
    difficult to do any funny business.

    Args:
        str_input (str): String to be evaluated.
        x (int|float): Fractional coordinate in :math:`x`.
        y (int|float): Fractional coordinate in :math:`y`.
        z (int|float): Fractional coordinate in :math:`z`.

    Returns:
        list[list[int|float,int|float,int|float]]:
            :math:`(N,3)` list of fractional coordinates.

    """
    ordered_inputs = {"x": "{0:.20f}", "y": "{1:.20f}", "z": "{2:.20f}"}
    # Replace any x, y, or z with the same character surrounded by curly braces. Then,
    # perform substitutions to insert the actual values.
    substituted_string = (
        re.sub(r"([xyz])", r"{\1}", str_input).format(**ordered_inputs).format(x, y, z)
    )

    # Remove any unexpected characters from the string.
    safe_string = re.sub(r"[^\d\[\]\,\+\-\/\*\.]", "", substituted_string)
    # Double check to be sure:
    assert all(
        char in ",.0123456789+-/*[]" for char in safe_string
    ), "Check that string only contains numerics or characters in { [],.+-/ }."
    return eval(safe_string, {"__builtins__": {}}, {})  # noqa: S307


def extract_unit_cell(filename: str, n_decimal_places: int = 4):
    """Return a complete unit cell from a CIF file in fractional coordinates.

    Args:
        filename (str): The name of the .cif file to be parsed.
        n_decimal_places (int, optional):
            The number of decimal places to round each position to for the uniqueness
            comparison. Because CIF files only store a few decimal places, a relatively
            low value is required for good results. 4 decimal places is usually enough
            to differentiate every unique position.
            Default value = ``4``

    Returns:
        :math:`(N, 3)` :class:`numpy.ndarray[np.float32]`:
            The full unit cell of the crystal structure.
    """
    fractional_positions = read_fractional_positions(filename=filename)

    symops = read_symmetry_operations(filename)
    symops_str = np.array2string(
        symops,
        separator=",",  # Place a comma after each line in the array. Required to eval
        threshold=np.inf,  # Ensure that every line is included in the string
        floatmode="unique",  # Ensures strings can uniquely represent each float number
    )

    all_frac_positions = [_safe_eval(symops_str, *xyz) for xyz in fractional_positions]

    pos = np.vstack(all_frac_positions)

    # Wrap particles into the box
    pos %= 1

    # TODO: add "fractional" flag?

    return np.unique(pos.round(n_decimal_places), axis=0)
