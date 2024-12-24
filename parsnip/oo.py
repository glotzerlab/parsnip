# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

r"""An interface for reading CIF files in Python.

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

        # A table with two columns and eight rows:
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

from __future__ import annotations

import re
import warnings
from collections.abc import Iterable

import numpy as np
from more_itertools import flatten, peekable
from numpy.lib.recfunctions import structured_to_unstructured

from parsnip._errors import ParseWarning
from parsnip.patterns import (
    _dtype_from_int,
    _is_data,
    _is_key,
    _line_is_continued,
    _safe_eval,
    _semicolon_to_string,
    _strip_comments,
    _strip_quotes,
    _try_cast_to_numeric,
    _write_debug_output,
    cast_array_to_float,
)
from parsnip.unitcells import _matrix_from_lengths_and_angles

# from parsnip.patterns import

NONTABLE_LINE_PREFIXES = ("_", "#")


class CifFile:
    """Parser for CIF files."""

    def __init__(self, fn: str, cast_values: bool = False):
        """Create a CifFile object from a filename.

        On construction, the entire file is parsed into key-value pairs and data tables.
        Comment lines are ignored.

        """
        self._fn = fn
        self._pairs = {}
        self._tables = []
        # self._table_keys = []
        # self._table_data = []

        self._cpat = {k: re.compile(pattern) for (k, pattern) in self.PATTERNS.items()}
        self._cast_values = cast_values

        with open(fn) as file:
            self._data = file.read()
        self._parse()

    @property
    def cast_values(self):
        """Whether to cast "number-like" values to ints & floats.

        .. note::

            When set to `True` after construction, the values are modified in-place.
            This action cannot be reversed.
        """
        return self._cast_values

    @cast_values.setter
    def cast_values(self, cast: bool):
        if cast:
            self._pairs = {
                k: _try_cast_to_numeric(_strip_quotes(v)) for (k, v) in self.pairs
            }
        else:
            warnings.warn(
                "Setting cast_values True->False has no effect on stored data.",
                category=ParseWarning,
                stacklevel=2,
            )
        self._cast_values = cast

    @property
    def pairs(self):
        """A dict containing key-value pairs extracted from the file.

        Numeric values will be parsed to int or float if possible. In these cases,
        precision specifiers will be stripped.

        Returns:
        --------
        dict[str : str | float | int]
        """
        return self._pairs

    @property
    def tables(self):
        """A list of data tables extracted from the file.

        These are stored as numpy structured arrays (see [the docs](https://numpy.org/doc/stable/user/basics.rec.html)
        for more information), which can be indexed by column labels.

        .. TODO helper function for converting to unstructured array


        Returns:
        --------
        list[(list[str], :math:`(N, N_{keys})` :class:`numpy.ndarray[str]`)]:
            A list of tuples corresponding with the keys and data of the system.
        """
        return self._tables

    @property
    def table_labels(self):
        """A list of column labels for each data array.

        This property is equivalent to `[arr.dtype.names for arr in self.tables]`.
        """
        return [arr.dtype.names for arr in self.tables]

    def _find_slice_in_tables(self, index: str):
        # TODO: only returns first match
        for table in self.tables:
            if index in table.dtype.names:
                return table[index]

    def get_from_tables(self, index: str | Iterable[str]):
        """Return a column or columns from the matching table in :meth:`~.self.tables`.

        If index is a single string, a single column will be returned from the matching
        table. If index is an Iterable of strings, the corresponding table slices will
        be returned. Slices from the same table will be grouped in the output array, but
        slices from different arrays will be returned seperately.

        .. tip::

            It is highly recommended to provide indices to seperate tables in seperate
            calls to this function.

            .. code:: python

                cf = CifFile(fn="my_file.cif")

                recommended_output = [
                    cf.get_from_tables(i) for i in [table_1_indices, table_2_indices]
                ]

                also_works = cf.get_from_tables([*table_1_indices, *table_2_indices])


        .. note::

            Returned arrays will match the ordering of input `index` keys *if all*
            *indices correspond to a single table.* Indices that match multiple tables
            will return all possible matches, in the order of the input tables. Lists of
            input that correspond with multiple tables will return data from those
            tables *in the order of the tablables in the :meth:`~.tables` property.*

        Parameters
        ----------
            index: str | Iterable[str]
                A column name or list of column names.

        Returns:
        --------
            list[:class:`numpy.ndarray`:] | :class:`numpy.ndarray`:
                A list of *unstructured* arrays corresponding with matches from the
                input keys. If the resulting list would have length 1, the data is
                returned directly instead. See the note above for data ordering.
        """
        index = np.atleast_1d(index)
        result = []
        for table in self.tables:
            matches = index[np.any(index[:, None] == table.dtype.names, axis=1)]
            if len(matches) == 0:
                continue

            result.append(
                structured_to_unstructured(
                    table[matches], copy=False, casting="safe"
                ).squeeze(axis=1)
            )
        return result if len(result) != 1 else result[0]

    PATTERNS = {
        "key_value_numeric": r"^(_[\w\.]+)[ |\t]+(-?\d+\.?\d*)",
        "key_value_general": r"^(_[\w\.-]+)[ |\t]+([^#^\n]+)",
        "table_delimiter": r"([Ll][Oo][Oo][Pp]_)[ |\t]*([^\n]*)",
        "block_delimiter": r"([Dd][Aa][Tt][Aa]_)[ |\t]*([^\n]*)",
        "key_list": r"_[\w_\.*]+[\[\d\]]*",
        "space_delimited_data": r"(\'[^\']*\'|\"[^\"]*\"]|[^\'\" \t]*)[ | \t]*",
    }

    def __getitem__(self, key: str | Iterable[str]):
        """Return an item from the dictionary of key-value pairs.

        Indexing with a string returns the value from the :meth:`~.pairs` dict. Indexing
        with an Iterable of strings returns a list of values, with `None` as a
        placeholder for keys that did not match any data.
        """
        if isinstance(key, Iterable) and not isinstance(key, str):
            return [self.pairs.get(k, None) for k in key]

        return self.pairs[key]

    def read_symmetry_operations(self):
        r"""Extract the symmetry operations from a CIF file.

        Args:
            filename (str): The name of the .cif file to be parsed.

        Returns:
            :math:`(N,)` :class:`numpy.ndarray[str]`:
                Symmetry operations as strings.
        """
        symmetry_keys = (
            "_symmetry_equiv_pos_as_xyz",
            "_space_group_symop_operation_xyz",
        )

        # Only one of the two keys will be matched. We can safely ignore that warning.
        warnings.filterwarnings("ignore", "Keys {'_", category=ParseWarning)
        return self.get_from_tables(symmetry_keys)

    def read_wyckoff_positions(self):
        r"""Extract symmetry-irreducible, fractional X,Y,Z coordinates from a CIF file.

        Parameters:
        -----------
            filename (str): The name of the .cif file to be parsed.

        Returns:
        --------
            :math:`(N, 3)` :class:`numpy.ndarray[np.float32]`:
                Fractional X,Y,Z coordinates of the unit cell.
        """
        xyz_keys = ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")
        xyz_data = self.get_from_tables(xyz_keys)
        xyz_data = cast_array_to_float(arr=xyz_data, dtype=np.float64)

        return xyz_data

    def read_cell_params(self, degrees: bool = True, mmcif: bool = False):
        r"""Read the cell lengths and angles from a CIF file.

        Paramters:
        ----------
            degrees (bool, optional):
                When True, angles are returned in degrees (as per the cif spec). When
                False, angles are converted to radians. Default value = ``True``
            mmcif (bool, optional):
                When False, the standard CIF key naming is used (e.g. _cell_angle_alpha)
                . When True, the mmCIF standard is used instead (e.g. cell.angle_alpha).
                Default value = ``False``

        Returns:
        --------
            tuple:
                The box vector lengths and angles in degrees or radians
                :math:`(L_1, L_2, L_3, \alpha, \beta, \gamma)`.
        """
        if mmcif:
            angle_keys = ("_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma")
            box_keys = (
                "_cell.length_a",
                "_cell.length_b",
                "_cell.length_c",
            ) + angle_keys
        else:
            angle_keys = ("_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")
            box_keys = (
                "_cell_length_a",
                "_cell_length_b",
                "_cell_length_c",
            ) + angle_keys
        cell_data = cast_array_to_float(arr=self[box_keys], dtype=np.float64)

        assert all(value is not None for value in cell_data)
        assert all(
            0 < key < 180 for key in cell_data[3:]
        ), "Read cell params were not in the expected range (0 < angle < 180 degrees)."

        if not degrees:
            cell_data[3:] = np.deg2rad(cell_data[3:])

        return tuple(cell_data)

    def extract_atomic_positions(
        self,
        fractional: bool = True,
        n_decimal_places: int = 4,
        verbose: bool = False,
    ):
        """Reconstruct atomic positions from Wyckoff sites and symmetry operations.

        .. warning::

            Reconstructing positions requires several floating point calculations that
            can be impacted by low-precision data in CIF files. Typically, at least four
            decimal places are required to accurately reconstruct complicated unit
            cells: less precision than this can yield cells with duplicate or missing
            positions.

        Args:
            fractional (bool, optional):
                Whether to return fractional or absolute coordinates.
                Default value = ``True``
            n_decimal_places (int, optional):
                The number of decimal places to round each position to for the
                uniqueness comparison. Values higher than 4 may not work for all CIF
                files. Default value = ``4``
            verbose (bool, optional):
                Whether to print debug information about the uniqueness checks.
                Default value = ``False``

        Returns:
            :math:`(N, 3)` :class:`numpy.ndarray[np.float32]`:
                The full unit cell of the crystal structure.
        """
        fractional_positions = self.read_wyckoff_positions()

        # Read the cell params and conver to a matrix of basis vectors
        cell = self.read_cell_params(degrees=False, mmcif=False)
        cell_matrix = _matrix_from_lengths_and_angles(*cell)

        symops = self.read_symmetry_operations()
        symops_str = np.array2string(
            symops,
            separator=",",  # Place a comma after each line in the array for eval
            threshold=np.inf,  # Ensure that every line is included in the string
            floatmode="unique",  # Ensures strings can uniquely represent each float
        )

        all_frac_positions = [
            _safe_eval(symops_str, *xyz) for xyz in fractional_positions
        ]

        pos = np.vstack(all_frac_positions)
        pos %= 1  # Wrap particles into the box

        # Filter unique points. This takes some time but makes the method faster overall
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

        return (
            pos[unique_indices] if fractional else real_space_positions[unique_indices]
        )

    def _parse(self):
        """Parse the cif file into python objects."""
        data_iter = peekable(self._data.split("\n"))

        for line in data_iter:
            if data_iter.peek(None) is None:
                break  # Exit without StopIteration

            # Combine nonsimple data values into a single, parseable line ==============
            while _line_is_continued(data_iter.peek()):
                line += _strip_comments(next(data_iter))
            line = _semicolon_to_string(line)

            # Skip processing if the line contains no data =============================
            if line == "" or _strip_comments(line) == "":
                continue

            # TODO: could support multi-block files in the future ======================
            block = re.match(self._cpat["block_delimiter"], line)
            if block is not None:
                continue

            # Extract key-value pairs and save to the internal state ===================
            pair = self._cpat["key_value_general"].match(line)
            if pair is not None:
                self._pairs.update(
                    {
                        pair.groups()[0]: _try_cast_to_numeric(
                            _strip_quotes(pair.groups()[1])
                        )
                        if self.cast_values
                        else pair.groups()[1].strip()
                    }
                )

            # Build up tables by incrementing through the iterator =====================
            table = re.match(self._cpat["table_delimiter"], line)

            if table is not None:
                table_keys, table_data = [], []

                # First, extract table headers. Must be prefixed with underscore
                line_groups = table.groups()
                if line_groups[-1] != "":  # Extract table keys from the _loop line
                    fragment = _strip_comments(line_groups[-1].strip())
                    if fragment[:1] == "_":
                        keys = self._cpat["key_list"].findall(fragment)
                        table_keys.extend([] if keys is None else keys)
                    else:
                        continue

                while _is_key(data_iter.peek(None)):
                    line = _strip_comments(next(data_iter))
                    while _line_is_continued(data_iter.peek(None)):
                        line += _strip_comments(next(data_iter))
                    line = _semicolon_to_string(line)
                    table_keys.extend(self._cpat["key_list"].findall(line))

                while _is_data(data_iter.peek(None)):
                    line = _strip_comments(next(data_iter))
                    while _line_is_continued(data_iter.peek(None)):
                        line += _strip_comments(next(data_iter))
                    line = _semicolon_to_string(line)
                    parsed_line = self._cpat["space_delimited_data"].findall(line)
                    parsed_line = [m for m in parsed_line if m != ""]
                    table_data.extend([parsed_line] if parsed_line else [])

                n_elements, n_cols = (
                    sum(len(row) for row in table_data),
                    len(table_keys),
                )

                if n_cols == 0:
                    continue  # Skip empty tables

                if n_elements % n_cols != 0:
                    warnings.warn(
                        f"Parsed data for table {len(self.tables)+1} cannot be resolved"
                        f" into a table of the expected size and will be ignored. "
                        f"Got n={n_elements} items, expected c={n_cols} columns: "
                        f"n%c={n_elements % n_cols}).",
                        category=ParseWarning,
                        stacklevel=2,
                    )
                    # self.tables.append((table_keys, table_data))
                    # self.tables.append(table_data)
                    continue
                if not all(len(key) == len(table_keys[0]) for key in table_keys):
                    table_data = np.array([*flatten(table_data)]).reshape(-1, n_cols)
                dt = _dtype_from_int(max(max(len(s) for s in l) for l in table_data))

                if len(set(table_keys)) < len(table_keys):
                    warnings.warn(
                        "Duplicate keys detected - table will not be processed.",
                        category=ParseWarning,
                        stacklevel=2,
                    )
                    continue

                rectable = np.atleast_2d(table_data)
                rectable.dtype = [*zip(table_keys, [dt] * n_cols)]
                rectable = rectable.reshape(rectable.shape, order="F")
                self.tables.append(rectable)

            if data_iter.peek(None) is None:
                break


if __name__ == "__main__":
    fn = "tests/sample_data/B-IncStrDb_Ccmm.cif"
    cf = CifFile(fn=fn)
    # [print(pair) for pair in cf.pairs.items()]
    # print(cf.tables[0])
