# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

r"""An interface for reading `CIF`_ files in Python.

.. _`CIF`: https://www.iucr.org/resources/cif

.. include:: ../../README.rst
    :start-after: .. _parse:
    :end-before: .. _installing:

.. admonition:: The CIF Format

    This is an example of a simple CIF file. A `key`_ (data name or tag) must start with
    an underscore, and is separated from the data value with whitespace characters.
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
from typing import ClassVar

import numpy as np
from more_itertools import flatten, peekable
from numpy.lib.recfunctions import structured_to_unstructured
from numpy.typing import ArrayLike

from parsnip._errors import ParseWarning
from parsnip.patterns import (
    _accumulate_nonsimple_data,
    _box_from_lengths_and_angles,
    _dtype_from_int,
    _is_data,
    _is_key,
    _matrix_from_lengths_and_angles,
    _safe_eval,
    _strip_comments,
    _strip_quotes,
    _try_cast_to_numeric,
    _write_debug_output,
    cast_array_to_float,
)

NONTABLE_LINE_PREFIXES = ("_", "#")


class CifFile:
    """Lightweight, performant parser for CIF files.

    Example
    -------
    To get started, simply provide a filename:

    >>> from parsnip import CifFile
    >>> cif = CifFile("example_file.cif")
    >>> print(cif)
    CifFile(fn=example_file.cif) : 12 data entries, 2 data loops

    Data entries are accessible via the :attr:`~.pairs` and :attr:`~.loops` attributes:

    >>> cif.pairs
    {'_journal_year': '1999', '_journal_page_first': '0', ...}
    >>> cif.loops[0]
    array([[('Cu1', '0.0000000000', '0.0000000000', '0.0000000000', 'Cu', 'a')]],
          dtype=...)
    >>> cif.loops[1]
    array([[('1', 'x,y,z')],
           [('96', 'z,y+1/2,x+1/2')],
           [('118', 'z+1/2,-y,x+1/2')],
           [('192', 'z+1/2,y+1/2,x')]],
          dtype=...)

    .. tip::

        See the docs for :attr:`__getitem__` and :attr:`get_from_loops` to query
        for data by key or column label!

    Parameters
    ----------
        fn : str
            Name of the file to be opened.
        cast_values : bool, optional
            Whether to convert string numerics to integers and float.
            Default value = ``False``
    """

    def __init__(self, fn: str, cast_values: bool = False):
        """Create a CifFile object from a filename.

        On construction, the entire file is parsed into key-value pairs and data loops.
        Comment lines are ignored.

        """
        self._fn = fn
        self._pairs = {}
        self._loops = []

        self._cpat = {k: re.compile(pattern) for (k, pattern) in self.PATTERNS.items()}
        self._cast_values = cast_values

        with open(fn) as file:
            self._parse(peekable(file))

    @property
    def pairs(self):
        """A dict containing key-value pairs extracted from the file.

        Numeric values will be parsed to int or float if possible. In these cases,
        precision specifiers will be stripped.

        Returns
        -------
        dict[str , str | float | int]
        """
        return self._pairs

    @property
    def loops(self):
        """A list of data tables (:code:``loop_``'s) extracted from the file.

        These are stored as `numpy structured arrays`_, which can be indexed by column
        labels. See the :attr:`~.structured_to_unstructured` helper function below for
        details on converting to standard arrays.

        .. _`numpy structured arrays`: https://numpy.org/doc/stable/user/basics.rec.html

        Returns
        -------
        list[:class:`numpy.ndarray[str]`]:
            A list of structured arrays containing table data from the file.
        """
        return self._loops

    @property
    def table_labels(self):
        """A list of column labels for each data array.

        This property is equivalent to :code:`[arr.dtype.names for arr in self.loops]`.

        Returns
        -------
        list[list[str]]:
            Column labels for :attr:`~.loops`, stored as a nested list of strings.
        """
        return [arr.dtype.names for arr in self.loops]

    def get_from_loops(self, index: ArrayLike):
        """Return a column or columns from the matching table in :attr:`~.loops`.

        If index is a single string, a single column will be returned from the matching
        table. If index is an Iterable of strings, the corresponding table slices will
        be returned. Slices from the same table will be grouped in the output array, but
        slices from different arrays will be returned separately.

        .. tip::

            It is highly recommended that queries across multiple loops are provided in
            separated calls to this function. This helps ensure output data is ordered
            as expected and allows for easier handling of cases where non-matching keys
            are provided.


        Example
        -------
        Extract a single column from a single table:

        >>> cif.get_from_loops("_symmetry_equiv_pos_as_xyz")
        array([['x,y,z'],
               ['z,y+1/2,x+1/2'],
               ['z+1/2,-y,x+1/2'],
               ['z+1/2,y+1/2,x']], dtype='<U14')

        Extract multiple columns from a single table:

        >>> table_1_cols = ["_symmetry_equiv_pos_site_id", "_symmetry_equiv_pos_as_xyz"]
        >>> cif.get_from_loops(table_1_cols)
        array([['1', 'x,y,z'],
               ['96', 'z,y+1/2,x+1/2'],
               ['118', 'z+1/2,-y,x+1/2'],
               ['192', 'z+1/2,y+1/2,x']], dtype='<U14')

        Extract multiple columns from multiple loops:

        >>> table_1_cols = ["_symmetry_equiv_pos_site_id", "_symmetry_equiv_pos_as_xyz"]
        >>> table_2_cols = ["_atom_site_type_symbol", "_atom_site_Wyckoff_label"]
        >>> [cif.get_from_loops(cols) for cols in (table_1_cols, table_2_cols)]
        [array([['1', 'x,y,z'],
               ['96', 'z,y+1/2,x+1/2'],
               ['118', 'z+1/2,-y,x+1/2'],
               ['192', 'z+1/2,y+1/2,x']], dtype='<U14'),
            array([['Cu', 'a']], dtype='<U12')]


        .. caution::

            Returned arrays will match the ordering of input ``index`` keys if all
            indices correspond to a single table. Indices that match multiple loops
            will return all possible matches, in the order of the input loops. Lists of
            input that correspond with multiple loops will return data from those
            loops *in the order they were read from the file.*

        Case where ordering of output matches the input file, not the provided keys:

        >>> cif.get_from_loops([*table_1_cols, *table_2_cols])
        [array([['Cu', 'a']], dtype='<U12'),
         array([['1', 'x,y,z'],
                ['96', 'z,y+1/2,x+1/2'],
                ['118', 'z+1/2,-y,x+1/2'],
                ['192', 'z+1/2,y+1/2,x']], dtype='<U14')]

        Parameters
        ----------
            index: str | Iterable[str]
                A column name or list of column names.

        Returns
        -------
            list[:class:`numpy.ndarray`] | :class:`numpy.ndarray`:
                A list of *unstructured* arrays corresponding with matches from the
                input keys. If the resulting list would have length 1, the data is
                returned directly instead. See the note above for data ordering.
        """
        index = np.atleast_1d(index)
        result = []
        for table in self.loops:
            matches = index[np.any(index[:, None] == table.dtype.names, axis=1)]
            if len(matches) == 0:
                continue

            result.append(
                structured_to_unstructured(
                    table[matches], copy=True, casting="safe"
                ).squeeze(axis=1)
            )
        return result if len(result) != 1 else result[0]

    def __getitem__(self, index: str | Iterable[str]):
        """Return an item from the dictionary of key-value pairs.

        Indexing with a string returns the value from the :meth:`~.pairs` dict. Indexing
        with an Iterable of strings returns a list of values, with `None` as a
        placeholder for keys that did not match any data.

        Example
        -------
        Indexing the class with a single key:

        >>> cif["_journal_year"]
        '1999'

        Indexing with a list of keys:

        >>> cif[["_journal_year", "_journal_page_first", "_journal_page_last"]]
        ['1999', '0', '123']

        Parameters
        ----------
            index: str | Iterable[str]
                An item key or list of keys.

        Returns
        -------
            list[str|int|float] :
                A list of data elements corresponding to the input key or keys. If the
                resulting list would have length 1, the item is returned directly
                instead.
        """
        if isinstance(index, Iterable) and not isinstance(index, str):
            return [self.pairs.get(k, None) for k in index]

        return self.pairs.get(index, None)

    def read_cell_params(self, degrees: bool = True, mmcif: bool = False):
        r"""Read the `unit cell parameters`_ (lengths and angles) from a CIF file.

        .. _`unit cell parameters`: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ccell.html

        Parameters
        ----------
            degrees : bool, optional
                When True, angles are returned in degrees (as per the CIF spec). When
                False, angles are converted to radians. Default value = ``True``
            mmcif : bool, optional
                When False, the standard CIF key naming is used (e.g. _cell_angle_alpha)
                . When True, the mmCIF standard is used instead (e.g. cell.angle_alpha).
                Default value = ``False``

        Returns
        -------
            tuple[float]:
                The box vector lengths in angstroms, and angles in degrees or radians
                :math:`(L_1, L_2, L_3, \alpha, \beta, \gamma)`.

        Raises
        ------
        ValueError
            If the stored data cannot form a valid box.
        """
        if mmcif:
            angle_keys = ("_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma")
            box_keys = (
                "_cell.length_a",
                "_cell.length_b",
                "_cell.length_c",
                *angle_keys,
            )
        else:
            angle_keys = ("_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")
            box_keys = (
                "_cell_length_a",
                "_cell_length_b",
                "_cell_length_c",
                *angle_keys,
            )
        if self.cast_values:
            cell_data = [float(x) for x in self[box_keys]]
        else:
            cell_data = cast_array_to_float(arr=self[box_keys], dtype=np.float64)

        def angle_is_invalid(x: float):
            return x <= 0.0 or x >= 180.0

        # if any(value is None for value in cell_data):
        #     missing = [k for k, v in zip(box_keys, cell_data) if v is None]
        #     msg = f"Keys {missing} did not return any data!"
        #     raise ValueError(msg) # TODO: reincorporate this error

        if any(angle_is_invalid(value) for value in cell_data[3:]):
            invalid = [
                k for k, v in zip(angle_keys, cell_data[3:]) if angle_is_invalid(v)
            ]
            msg = f"Keys {invalid} are outside the valid range (0 <= θ <= 180)."
            raise ValueError(msg)

        if not degrees:
            cell_data[3:] = np.deg2rad(cell_data[3:])

        return tuple(float(v) for v in cell_data)  # Return as base python types

    def read_symmetry_operations(self):
        r"""Extract the symmetry operations from a CIF file.

        Returns
        -------
            :math:`(N,)` :class:`numpy.ndarray[str]`:
                An array of strings containing the symmetry operations in a
                `parsable algebraic form`_.

        .. _`parsable algebraic form`: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_symop_operation_xyz.html
        """
        symmetry_keys = (
            "_symmetry_equiv_pos_as_xyz",
            "_space_group_symop_operation_xyz",
        )

        # Only one key is valid in each standard, so we only ever get one match.
        return self.get_from_loops(symmetry_keys)

    def read_wyckoff_positions(self):
        r"""Extract symmetry-irreducible, fractional x,y,z coordinates from a CIF file.

        Returns
        -------
            :math:`(N, 3)` :class:`numpy.ndarray[float]`:
                Symmetry-irreducible positions of atoms in `fractional coordinates`_.

        .. _`fractional coordinates`: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html
        """
        xyz_keys = ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")

        return cast_array_to_float(arr=self.get_from_loops(xyz_keys), dtype=float)

    def build_unit_cell(
        self,
        fractional: bool = True,
        n_decimal_places: int = 4,
        wrap_coords: bool = True,  # TODO: docs
        verbose: bool = False,
    ):
        """Reconstruct atomic positions from Wyckoff sites and symmetry operations.

        Rather than storing an entire unit cell's atomic positions, CIF files instead
        include the data required to recreate those positions based on symmetry rules.
        Symmetry operations (stored as strings of x,y,z position permutations) are
        applied to the Wyckoff (symmetry irreducible) positions to create a list of
        possible atomic sites. These are then wrapped into the unit cell and filtered
        for uniqueness to yield the final crystal.

        .. warning::

            Reconstructing positions requires several floating point calculations that
            can be impacted by low-precision data in CIF files. Typically, at least four
            decimal places are required to accurately reconstruct complicated unit
            cells: less precision than this can yield cells with duplicate or missing
            positions.

        Args:
            fractional : (bool, optional)
                Whether to return fractional or absolute coordinates.
                Default value = ``True``
            n_decimal_places : (int, optional)
                The number of decimal places to round each position to for the
                uniqueness comparison. Values higher than 4 may not work for all CIF
                files. Default value = ``4``
            verbose : (bool, optional)
                Whether to print debug information about the uniqueness checks.
                Default value = ``False``

        Returns
        -------
            :math:`(N, 3)` :class:`numpy.ndarray[float]`:
                The full unit cell of the crystal structure.

        Raises
        ------
        ValueError
            If the stored data cannot form a valid box.
        """
        fractional_positions = self.read_wyckoff_positions()

        # Read the cell params and convert to a matrix of basis vectors
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
        if wrap_coords:
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

    @property
    def cast_values(self):
        """Bool : Whether to cast "number-like" values to ints & floats.

        .. note::

            When set to `True` after construction, the values are modified in-place.
            This action cannot be reversed.
        """
        return self._cast_values

    @cast_values.setter
    def cast_values(self, cast: bool):
        if cast:
            self._pairs = {
                k: _try_cast_to_numeric(_strip_quotes(v))
                for (k, v) in self.pairs.items()
            }
        else:
            warnings.warn(
                "Setting cast_values True->False has no effect on stored data.",
                category=ParseWarning,
                stacklevel=2,
            )
        self._cast_values = cast

    @property
    def cell(self):
        """Read the unit cell as a `freud`_ or HOOMD `box-like`_ object.

        .. _`box-like`: https://hoomd-blue.readthedocs.io/en/v5.0.0/hoomd/module-box.html#hoomd.box.box_like
        .. _`freud`: https://freud.readthedocs.io/en/latest/gettingstarted/examples/module_intros/box.Box.html

        Example
        -------
        This method provides a convinient interface to create box objects.

        >>> import freud
        >>> cell = cif.cell
        >>> print(cell)
        (3.6, 3.6, 3.6, 0.0, 0.0, 0.0)
        >>> freud.Box(*cell) # doctest: +SKIP
        freud.box.Box(Lx=3.6, Ly=3.6, Lz=3.6, xy=0, xz=0, yz=0, ...)
        >>> hoomd.Box(*cell) # doctest: +SKIP
        hoomd.box.Box(Lx=3.6, Ly=3.6, Lz=3.6, xy=0.0, xz=0.0, yz=0.0)


        Returns
        -------
        tuple:
            The box vector lengths (in Ångströms) and tilt factors (unitless).
            :math:`(L_1, L_2, L_3, xy, xz, yz)`.
        """
        return _box_from_lengths_and_angles(
            *self.read_cell_params(degrees=False, mmcif=False)
        )

    @classmethod
    def structured_to_unstructured(cls, arr: np.ndarray):
        """Convert a structured (column-labeled) array to a standard unstructured array.

        This is useful when extracting entire loops from :attr:`~.loops` for use in
        other programs. This classmethod simply calls
        :code:`np.lib.recfunctions.structured_to_unstructured` on the input data to
        ensure the resulting array is properly laid out in memory. See
        `this page in the structured array docs`_ for more information.

        .. _`this page in the structured array docs`: https://numpy.org/doc/stable/user/basics.rec.html

        Parameters
        ----------
            arr : :class:`numpy.ndarray`: | :class:`numpy.recarray`
                The structured array to convert.

        Returns
        -------
            :class:`numpy.ndarray`:
                An *unstructured* array containing a copy of the data from the input.
        """
        return structured_to_unstructured(arr, copy=True)

    def _parse(self, data_iter: Iterable):
        """Parse the cif file into python objects."""
        for line in data_iter:
            if data_iter.peek(None) is None:
                break  # Exit without StopIteration

            # Combine nonsimple data entries into a single, parseable line =============
            line = _accumulate_nonsimple_data(data_iter, line)

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
                        else pair.groups()[1].rstrip()  # Skip trailing newlines
                    }
                )

            # Build up tables by incrementing through the iterator =====================
            loop = re.match(self._cpat["loop_delimiter"], line)

            if loop is not None:
                loop_keys, loop_data = [], []

                # First, extract table headers. Must be prefixed with underscore
                line_groups = loop.groups()
                if line_groups[-1] != "":  # Extract loop keys from the _loop line
                    fragment = _strip_comments(line_groups[-1].strip())
                    if fragment[:1] == "_":
                        keys = self._cpat["key_list"].findall(fragment)
                        loop_keys.extend(keys if keys is not None else [])
                    else:
                        continue

                while _is_key(data_iter.peek(None)):
                    line = _accumulate_nonsimple_data(
                        data_iter, _strip_comments(next(data_iter))
                    )
                    loop_keys.extend(self._cpat["key_list"].findall(line))

                while _is_data(data_iter.peek(None)):
                    line = _accumulate_nonsimple_data(
                        data_iter, _strip_comments(next(data_iter))
                    )
                    parsed_line = self._cpat["space_delimited_data"].findall(line)
                    parsed_line = [m for m in parsed_line if m != ""]
                    loop_data.extend([parsed_line] if parsed_line else [])

                n_elements, n_cols = (
                    sum(len(row) for row in loop_data),
                    len(loop_keys),
                )

                if n_cols == 0:
                    continue  # Skip empty tables

                if n_elements % n_cols != 0:
                    warnings.warn(
                        f"Parsed data for table {len(self.loops)+1} cannot be resolved"
                        f" into a table of the expected size and will be ignored. "
                        f"Got n={n_elements} items, expected c={n_cols} columns: "
                        f"n%c={n_elements % n_cols}).",
                        category=ParseWarning,
                        stacklevel=2,
                    )
                    continue
                if not all(len(key) == len(loop_keys[0]) for key in loop_keys):
                    loop_data = np.array([*flatten(loop_data)]).reshape(-1, n_cols)
                dt = _dtype_from_int(max(max(len(s) for s in l) for l in loop_data))

                if len(set(loop_keys)) < len(loop_keys):
                    warnings.warn(
                        "Duplicate keys detected - table will not be processed.",
                        category=ParseWarning,
                        stacklevel=2,
                    )
                    continue

                rectable = np.atleast_2d(loop_data)
                rectable.dtype = [*zip(loop_keys, [dt] * n_cols)]
                rectable = rectable.reshape(rectable.shape, order="F")
                self.loops.append(rectable)

            if data_iter.peek(None) is None:
                break

    def __repr__(self):
        n_pairs = len(self.pairs)
        n_tabs = len(self.loops)
        return f"CifFile(fn={self._fn}) : {n_pairs} data entries, {n_tabs} data loops"

    PATTERNS: ClassVar = {
        "key_value_general": r"^(_[\w\.\-/\[\d\]]+)\s+([^#]+)",
        "loop_delimiter": r"([Ll][Oo][Oo][Pp]_)[ |\t]*([^\n]*)",
        "block_delimiter": r"([Dd][Aa][Tt][Aa]_)[ |\t]*([^\n]*)",
        "key_list": r"_[\w_\.*]+[\[\d\]]*",
        "space_delimited_data": r"(\;[^\;]*\;|\'[^\']*\'|\"[^\"]*\"]|[^\'\"\;\s]*)\s*",
    }
    """Regex patterns used when parsing files.

    This dictionary can be modified to change parsing behavior, although doing is not
    recommended. Changes to this variable are shared across all instances of the class.
    """
