"""Functions for constructing unit cells from CIF file data.

Rather than storing an entire unit cell's atomic positions, CIF files instead include
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


def _matrix_from_lengths_and_angles(l1, l2, l3, alpha, beta, gamma):
    a1 = np.array([l1, 0, 0])
    a2 = np.array([l2 * np.cos(gamma), l2 * np.sin(gamma), 0])
    a3x = np.cos(beta)
    a3y = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    under_sqrt = 1 - a3x**2 - a3y**2
    if under_sqrt < 0:
        raise ValueError("The provided angles can not form a valid box.")
    a3z = np.sqrt(under_sqrt)
    a3 = np.array([l3 * a3x, l3 * a3y, l3 * a3z])
    return np.array([a1, a2, a3]).T


def read_wyckoff_positions(
    filename: str,
    regex_filter: tuple[tuple[str, str]] | None = ((r",\s+", ",")),
):
    r"""Extract the symmetry-irreducible, fractional X,Y,Z coordinates from a CIF file.

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


def _write_debug_output(unique_indices, unique_counts, pos, check="Initial"):
    print(f"{check} uniqueness check:")
    if len(unique_indices) == len(pos):
        print("... all points are unique (within tolerance).")
    else:
        print("(duplicate point, number of occurences)")
        [
            print(pt, count)
            for pt, count in zip(pos[unique_indices], unique_counts)
            if count > 1
        ]

    print()


def extract_atomic_positions(
    filename: str,
    fractional: bool = True,
    n_decimal_places: int = 4,
    verbose: bool = False,
):
    """Reconstruct atomic positions from Wyckoff sites and symmetry operations.

    .. warning::

        Reconstructing positions requires several floating point calculations that can
        be impacted by low-precision data in CIF files. Typically, at least four decimal
        places are required to accurately reconstruct complicated unit cells: less
        precision than this can yield cells with duplicate or missing positions.

    Args:
        filename (str): The name of the .cif file to be parsed.
        fractional (bool, optional):
            Whether to return fractional or absolute coordinates.
            Default value = ``True``
        n_decimal_places (int, optional):
            The number of decimal places to round each position to for the uniqueness
            comparison. Values higher than 4 may not work for all CIF files.
            Default value = ``4``
        verbose (bool, optional):
            Whether to print debug information about the uniqueness checks.
            Default value = ``False``

    Returns:
        :math:`(N, 3)` :class:`numpy.ndarray[np.float32]`:
            The full unit cell of the crystal structure.
    """
    fractional_positions = read_wyckoff_positions(filename=filename)

    # Read the cell params and conver to a matrix of basis vectors
    cell = read_cell_params(filename, degrees=False, mmcif=False)
    cell_matrix = _matrix_from_lengths_and_angles(*cell)

    symops = read_symmetry_operations(filename)
    symops_str = np.array2string(
        symops,
        separator=",",  # Place a comma after each line in the array. Required for eval
        threshold=np.inf,  # Ensure that every line is included in the string
        floatmode="unique",  # Ensures strings can uniquely represent each float number
    )

    all_frac_positions = [_safe_eval(symops_str, *xyz) for xyz in fractional_positions]

    pos = np.vstack(all_frac_positions)
    pos %= 1  # Wrap particles into the box

    # Filter unique points. This takese some time, but makes the method faster overall
    _, unique_indices, unique_counts = np.unique(
        pos.round(n_decimal_places), return_index=True, return_counts=True, axis=0
    )

    if verbose:
        _write_debug_output(unique_indices, unique_counts, pos, check="Initial")

    # Remove initial duplicates, then map to real space for a second check
    pos = pos[unique_indices]
    real_space_positions = pos @ cell_matrix

    _, unique_indices, unique_counts = np.unique(
        real_space_positions.round(n_decimal_places),
        return_index=True,
        return_counts=True,
        axis=0,
    )

    if verbose:
        _write_debug_output(unique_indices, unique_counts, pos, check="Secondary")

    """
    # This code allows for parity with Gemmi - however, the results are effectively
    # identical to the code above. This could be re-enabled in the future if desired.
    if not np.isclose(merge_dist,0):
        dists,i_inds,j_inds = _get_distances(real_space_positions[unique_indices])

        # Now, get the positions that are less than merge_dist apart and remove them
        overlapping_point_indices = np.vstack([i_inds,j_inds])[:, dists<merge_dist**2].T
        pos = np.delete(pos, overlapping_point_indices[:,1],axis=0)
        if verbose:
            print("Tertiary uniqueness check:")
            print(f"... {len(overlapping_point_indices)} points removed")
        real_space_positions = pos@cell_matrix
    """

    return pos[unique_indices] if fractional else real_space_positions[unique_indices]
