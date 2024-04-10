r"""Functions for parsing CIF files in Python.

.. include:: ../../README.rst
    :start-after: .. _parse:
    :end-before: .. _installing:

.. admonition:: The CIF Format

    This is an example of a simple CIF file. A `key`_ (data name or tag) must start with
    an underscore, and is seperated from the data value with whitespace characters.
    A `table`_ begins with the ``loop_`` keyword, and contain a header block and a data
    block. The vertical position of a tag in the table headings corresponds with the
    horizontal position of the associated column in the table values.

    .. code-block:: text

        # Key-value pairs describing the unit cell:
        _cell_length_a  5.40
        _cell_length_b  3.43
        _cell_length_c  5.08
        _cell_angle_alpha  90.0
        _cell_angle_beta  132.3
        _cell_angle_gamma  90.0

        # A table with two columns and three rows:
        loop_
        _symmetry_equiv_pos_site_id
        _symmetry_equiv_pos_as_xyz
        1  x,y,z
        2  -x,y,-z
        3  -x,-y,-z
        4  x,-y,z
        5  x+1/2,y+1/2,z
        6  -x+1/2,y+1/2,-z
        7  -x+1/2,-y+1/2,-z
        8  x+1/2,-y+1/2,z

        _symmetry_space_group_name_H-M  'C2 / m' # One more key-value pair


.. _key: https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#definitions
.. _table: https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#onelevel

"""
# TODO: Move "special" readers to seperate module?

from __future__ import annotations

import re
import warnings

import numpy as np

from ._errors import ParseError, ParseWarning
from ._utils import _deg2rad, _str2num
from .patterns import LineCleaner, cast_array_to_float, remove_nondelimiting_whitespace


def _remove_comments_from_line(line):
    return line.split("#")[0].strip()


def read_table(
    filename: str,
    keys: str,
    nondelimiting_whitespace_replacement: str = "_",
    regex_filter: tuple[str, str] | None = None,
    keep_original_key_order=False,
) -> np.ndarray[str]:
    r"""Extract data from a CIF file loop\_ table.

    .. tip::

        CIF tables are whitespace delimited - however, values enclosed in quotation
        marks may also contain whitespace characters. The parameter
        ``nondelimiting_whitespace_replacement`` handles this possibility by replacing
        nondelimiting whitespaces with underscores. This value can be also be set to an
        empty string, or any arbitrary sequence of characters.

    Args:
        filename (str):
            The name of the .cif file to be parsed.
        keys (tuple[str]):
            The names of the keys to be parsed. The columns associated with these keys
            will be returned in the final array.
        nondelimiting_whitespace_replacement (str, optional):
            Character to replace non-delimiting whitespaces with. By default,
        regex_filter (tuple[str,str], optional):
            A tuple of strings that are compiled to a regex filter and applied to each
            data line. If a tuple of tuples of strings is provided instead, each pattern
            will be applied seperately.
            Default value = ``None``
        keep_original_key_order (bool, optional):
            When True, preserve the order of keys in the table from the cif file.
            When False, return columns of data in order of the input ``keys`` arg.
            Default value = ``False``


    Returns:
        :math:`(N, N_{keys})` :class:`numpy.ndarray[str]`:
            A numpy array of the data as strings.

    .. warning::

        This function will ONLY return data from a single table. If keys are provided
        that correspond to data from multiple tables, only the first table will be read.

    .. tip::
        The ``filter_line`` argument allows for dynamic input creation of regex filters
        to apply to each line that contains data to be saved. The first value in the
        tuple is the pattern to match, and the second value is the replacement text.
        The default value is ``((",\s+",","))``, which removes whitespace following
        commas to help differentiate between individual data entries that contain spaces
        from other sections of the line that are also whitespace separated. Adding
        another tuple to remove single quotes can also be helpful: try
        ``((",\s+",","),("'",""))`` to achieve this. To disable the feature entirely,
        pass in a tuple of empty strings: ``("","")``.

    """
    with open(filename) as f:
        tables = f.read().split("loop_")

    if regex_filter is not None:
        line_cleaner = LineCleaner(regex_filter)
    nontable_line_prefixes = ("_", "#")

    for table in tables:
        lines = table.strip().split("\n")
        in_header = True
        data_column_indices, data, column_order = [], [], []

        for line_number, line in enumerate(lines):
            # Check for invalid blank lines in the table header
            if in_header and data_column_indices and line == "":
                raise ParseError(
                    "Whitespace may not be used in between keys in the table header. "
                    "See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#general"
                    ", section 7 for more details."
                )

            # We will get errors if there is a comment after the loop_ block that
            # contains our data. This is questionably legal, but very uncommon

            line = _remove_comments_from_line(line)

            # Save current key position if it is one of the keys we want.
            if in_header and (line in keys):
                data_column_indices.append(line_number)
                if not keep_original_key_order:
                    column_order.append(keys.index(line))
                continue

            # If we exit the header and enter the table body
            if data_column_indices and (line[:1] not in nontable_line_prefixes):
                in_header = False  # Exit the header and start writing data

                if regex_filter is not None:  # Apply user-defined regex, if present
                    line = line_cleaner(line)

                clean_line = remove_nondelimiting_whitespace(
                    line.strip(), replacement=nondelimiting_whitespace_replacement
                )
                split_line = clean_line.split()

                # Only add data if the line has at least as many columns as required.
                n_cols_found, n_cols_expected = (
                    len(split_line),
                    len(data_column_indices),
                )
                if n_cols_found >= n_cols_expected:
                    data.append(split_line)
                elif split_line != [] and n_cols_found < n_cols_expected:
                    warnings.warn(
                        f"Data line is a fragment and will be skipped: (expected line "
                        f"with {n_cols_expected} values, got {split_line}).",
                        ParseWarning,
                        stacklevel=2,
                    )
                continue
            elif (not in_header) and (line[:1] == "_"):
                break
        if data_column_indices:
            break

    if not keep_original_key_order:
        # Reorder the column indices to match the order of the input keys
        data_column_indices = np.array(data_column_indices)[np.argsort(column_order)]

    if len(column_order) != len(keys):
        missing_keys = {key for i, key in enumerate(keys) if i not in column_order}
        warnings.warn(
            f"Keys {missing_keys} were not found in the table.",
            ParseWarning,
            stacklevel=2,
        )
    return np.atleast_2d(data)[:, data_column_indices]


def _parsed_line_generator(filename, regexp):
    """Apply a regex pattern line by line and yield the pattern's matches.

    This is intended to be an internal function that handles the reading of CIF files.
    Abstracting this out clarifies which logic belongs to the file parser and which
    belongs to the actual data manipulation.

    Args:
        filename (str): The name of the .cif file to be parsed.
        regexp (str): tring to generate the regex pattern that is applied to each line.

    Yields:
        tuple(str,str|float|int):
    """
    pattern = re.compile(regexp)
    with open(filename) as file:
        for line in file:
            # Line is either empty, or does not start with a valid key
            if line == "" or line[0] != "_":
                continue
            parsed_line = pattern.match(line)
            if parsed_line:  # Regex matches
                yield parsed_line


def read_key_value_pairs(
    filename: str,
    keys: tuple[str] | None = None,
    only_read_numerics: bool = False,
):
    """Extract key-value pairs from a CIF file.

    By default, this function reads all keys and returns data values as strings. Setting
    ``only_read_numerics`` to True will cause the program to cast data to a nunmeric
    type (float or int): however, keys that cannot be safely cast into a numeric are
    skipped.

    Args:
        filename (str): The name of the .cif file to be parsed.
        keys (tuples[str]|None, optional):
            A tuple of keys to search and return data for.
            If keys is None, all keys are returned (Default value: None).
        only_read_numerics (bool, optional):
            Whether to read only values that cannot be cast to int or float.
            (Default value: False)

    Returns:
        dict[str,float|int] | dict[str,str]:
            Dict of the key value pairs. Values will either be all strings, or a mixture
            of int and float, and the order will match the order of keys (if provided).

    .. note::

            If no data is found for any of the provided keys, a warning will be raised
            and the output value will by ``None``.

    """
    # REGEX EXPLANATION
    # ^         : Match only at the start of the line
    # (_[\w-]+) : Match any number/mix of alphanumerics, "-", and "_", as a group
    # [ |\t]+   : Match one or more whitespace " " or tab characters.

    # Parse numbers:
    # (         : Start new group
    # -?\d+     : Match 0 or 1 "-" characters, then 1 or more digits 0-9
    # \.?       : Match 0 or 1 "." characters
    # \d*       : Match 0 or more digits 0-9
    # )         : End the group

    # Parse strings:
    # (         : Start new group
    # [^#^\n]+  : Match 1 ore more characters that are NOT a "#" or newline "\n"
    # )         : End the group

    data = {}

    if only_read_numerics:
        regexp = r"^(_\w+)[ |\t]+(-?\d+\.?\d*)"
    else:
        regexp = r"^(_[\w-]+)[ |\t]+([^#^\n]+)"

    if keys is not None:
        # Insertion order our dict with original key order
        for key in keys:
            data[key] = None
        # Convert to mutable datastructure so we can remove identified keys
        keys = set(keys)

    for parsed_line in _parsed_line_generator(filename, regexp=regexp):
        key, val = parsed_line.groups()
        val = _str2num(val) if only_read_numerics else val.strip()

        if keys is None:
            data[key] = val
        elif key in keys:
            data[key] = val
            keys.remove(key)
        elif len(keys) == 0:
            break

    if keys is not None and len(keys) != 0:
        warnings.warn(
            f"Keys {keys} did not match any data!", ParseWarning, stacklevel=2
        )

    return data


def read_cell_params(filename, degrees=True, validate=True):
    r"""Read the cell lengths and angles from a CIF file.

    Args:
        filename (str): The name of the .cif file to be parsed.
        degrees (bool, optional):
            When True, angles are returned in degrees (as per the cif spec). When False,
            angles are converted to radians (Default value: True).
        validate (bool, optional):
            Whether to check if the results are correct. (Default value: True)

    Returns:
        tuple:
            The box vector lengths and angles in degrees or radians
            :math:`(L_1, L_2, L_3, \alpha, \beta, \gamma)`.
    """
    angle_keys = ("_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")
    box_keys = ("_cell_length_a", "_cell_length_b", "_cell_length_c") + angle_keys

    cell_data = read_key_value_pairs(filename, keys=box_keys, only_read_numerics=True)

    if validate:
        assert all(value is not None for value in cell_data.values())
        assert all(0 < cell_data[key] < 180 for key in angle_keys)

    if not degrees:
        for key in angle_keys:
            cell_data[key] = _deg2rad(cell_data[key])

    return tuple(cell_data.values())


def read_fractional_positions(
    filename: str,
    filter_line: tuple[tuple[str, str]] = ((r",\s+", ",")),
):
    r"""Extract the fractional X,Y,Z coordinates from a CIF file.

    Args:
        filename (str): The name of the .cif file to be parsed.
        filter_line (tuple[tuple[str]], optional):
            A tuple of strings that are compiled to a regex filter and applied to each
            data line. (Default value: ((r",\s+",",")) )

    Returns:
        :math:`(N, 3)` :class:`numpy.ndarray[np.float32]`:
            Fractional X,Y,Z coordinates of the unit cell.
    """
    xyz_keys = ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")
    # Once #6 is added, we should warnings.catch_warnings(action="error")
    xyz_data = read_table(
        filename=filename,
        keys=xyz_keys,
    )

    xyz_data = cast_array_to_float(arr=xyz_data, dtype=np.float32)

    # Validate results
    assert xyz_data.shape[1] == 3
    assert xyz_data.dtype == np.float32

    return xyz_data
