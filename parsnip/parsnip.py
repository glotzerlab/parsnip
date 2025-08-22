# Copyright (c) 2025, The Regents of the University of Michigan
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

        # A header describing this portion of the file
        data_cif_Cu-FCC

        # Several key-value pairs
        _journal_year 1999
        _journal_page_first 0
        _journal_page_last 123

        _chemical_name_mineral 'Copper FCC'
        _chemical_formula_sum 'Cu'

        # Key-value pairs describing the unit cell (Å and °)
        _cell_length_a     3.6
        _cell_length_b     3.6
        _cell_length_c     3.6
        _cell_angle_alpha  90.0
        _cell_angle_beta   90.0
        _cell_angle_gamma  90.0

        # A table with 6 columns and one row
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        _atom_site_type_symbol
        _atom_site_Wyckoff_label
        Cu1 0.0000000000 0.0000000000 0.0000000000  Cu a

        _symmetry_space_group_name_H-M  'Fm-3m' # One more key-value pair

        # A table with two columns and four rows:
        loop_
        _symmetry_equiv_pos_site_id
        _symmetry_equiv_pos_as_xyz
        1  x,y,z
        96  z,y+1/2,x+1/2
        118  z+1/2,-y,x+1/2
        192  z+1/2,y+1/2,x


.. _key: https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#definitions
.. _table: https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#onelevel

"""

from __future__ import annotations

import re
import warnings
from collections.abc import Iterable
from fnmatch import filter as fnfilter
from fnmatch import fnmatch
from importlib.util import find_spec
from pathlib import Path
from typing import ClassVar

import numpy as np
from more_itertools import flatten, peekable
from numpy.lib.recfunctions import structured_to_unstructured

from parsnip._errors import ParseWarning
from parsnip.patterns import (
    _CIF_KEY,
    _PROG_PLUS,
    _PROG_STAR,
    _WHITESPACE,
    _accumulate_nonsimple_data,
    _box_from_lengths_and_angles,
    _bracket_pattern,
    _dtype_from_int,
    _flatten_or_none,
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
        fn : str | Path
            Path to the file to be opened.
        cast_values : bool, optional
            Whether to convert string numerics to integers and float.
            Default value = ``False``
    """

    def __init__(self, fn: str | Path, cast_values: bool = False):
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

    _SYMPY_AVAILABLE = find_spec("sympy") is not None

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
        r"""A list of data tables (`loop_`'s) extracted from the file.

        These are stored as `numpy structured arrays`_, which can be indexed by column
        labels. See the :attr:`~.structured_to_unstructured` helper function below for
        details on converting to standard arrays.

        .. _`numpy structured arrays`: https://numpy.org/doc/stable/user/basics.rec.html

        Returns
        -------
        list[numpy.ndarray[str]]
            A list of structured arrays containing table data from the file.
        """
        return self._loops

    def __getitem__(self, index: str | Iterable[str]):
        """Return an item or list of items from :meth:`~.pairs` and :meth:`~.loops`.

        This getter searches the entire CIF state to identify the input keys, returning
        ``None`` if the key does not match any data. Matching columns from `loop` tables
        are returned as 1D arrays.

        .. tip::

            This method of accessing data is recommended for most uses, as it ensures
            data is returned wherever possible. :meth:`~.get_from_loops` may be useful
            when multi-column slices of an array are needed.

        Example
        -------
        Indexing the class with a single key:

        >>> cif["_journal_year"]
        '1999'
        >>> cif["_atom_site_label"]
        array([['Cu1']], dtype='<U12')

        Indexing with a list of keys:

        >>> cif[["_chemical_name_mineral", "_symmetry_equiv_pos_as_xyz"]]
        ["'Copper FCC'",
        array([['x,y,z'],
            ['z,y+1/2,x+1/2'],
            ['z+1/2,-y,x+1/2'],
            ['z+1/2,y+1/2,x']], dtype='<U14')]

        Wildcards are supported for lookups with this method:

        >>> cif[["_journal*", "_atom_site_fract_?"]]
        [['1999', '0', '123'],
        ...array([['0.0000000000', '0.0000000000', '0.0000000000']], dtype='<U12')]

        Parameters
        ----------
            index: str | typing.Iterable[str]
                An item key or list of keys.
        """
        output = []
        index = [index] if isinstance(index, str) else index
        for key in index:
            pairs_match = self.get_from_pairs(key)
            loops_match = self.get_from_loops(key)
            output.append(pairs_match if pairs_match is not None else loops_match)
        return output[0] if len(output) == 1 else output

    def get_from_pairs(self, index: str | Iterable[str]):
        """Return an item or items from the dictionary of key-value pairs.

        .. tip::

            This method supports a few unix-style wildcards. Use ``*`` to match any
            number of any character, and ``?`` to match any single character. If a
            wildcard matches more than one key, a list is returned for that index.
            Lookups using this method are case-insensitive, per the CIF spec.

        Indexing with a string returns the value from the :meth:`~.pairs` dict. Indexing
        with an Iterable of strings returns a list of values, with ``None`` as a
        placeholder for keys that did not match any data.

        Example
        -------
        Indexing the class with a single key:

        >>> cif.get_from_pairs("_journal_year")
        '1999'

        Indexing with a list of keys:

        >>> cif.get_from_pairs(["_journal_page_first", "_journal_page_last"])
        ['0', '123']

        Indexing with wildcards:

        >>> cif.get_from_pairs("_journal*")
        ['1999', '0', '123']

        Single-character wildcards can generalize keys across CIF and mmCIF files:

        >>> cif.get_from_pairs("_symmetry?space_group_name_H-M")
        "'Fm-3m'"

        Parameters
        ----------
            index: str | typing.Iterable[str]
                An item key or list of keys.

        Returns
        -------
            list[str|int|float] :
                A list of data elements corresponding to the input key or keys. If the
                resulting list would have length 1, the item is returned directly
                instead.
        """
        if isinstance(index, str):  # Escape brackets with []
            index = re.sub(_bracket_pattern, r"[\1]", index)
            return _flatten_or_none(
                [
                    v
                    for (k, v) in self.pairs.items()
                    if fnmatch(k.lower(), index.lower())
                ]
            )

        # Escape all brackets in all indices
        index = [re.sub(_bracket_pattern, r"[\1]", i) for i in index]
        matches = [
            [
                _flatten_or_none(
                    [
                        v
                        for (k, v) in self.pairs.items()
                        if fnmatch(k.lower(), pat.lower())
                    ]
                )
            ]
            for pat in index
        ]
        return [_flatten_or_none(m) for m in matches]

    def get_from_loops(self, index: str | Iterable[str]):
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
        >>> table_1 = cif.get_from_loops(table_1_cols)
        >>> table_1
        array([['1', 'x,y,z'],
               ['96', 'z,y+1/2,x+1/2'],
               ['118', 'z+1/2,-y,x+1/2'],
               ['192', 'z+1/2,y+1/2,x']], dtype='<U14')

        Wildcard patterns are accepted for single input keys:

        >>> assert (cif.get_from_loops("_symmetry_equiv_pos*") == table_1).all()

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
            index: str | typing.Iterable[str]
                A column name or list of column names.

        Returns
        -------
            list[:class:`numpy.ndarray`] | :class:`numpy.ndarray`:
                A list of *unstructured* arrays corresponding with matches from the
                input keys. If the resulting list would have length 1, the data is
                returned directly instead. See the note above for data ordering.
        """
        if isinstance(index, str):
            result, index = [], re.sub(_bracket_pattern, r"[\1]", index)
            for table, labels in zip(self.loops, self.loop_labels):
                match = table[fnfilter(labels, index)]
                if match.size > 0:
                    result.append(
                        structured_to_unstructured(
                            match, copy=True, casting="safe"
                        ).squeeze(axis=1)
                    )
            if result == [] or (len(result) == 1 and len(result[0]) == 0):
                return None
            return result[0] if len(result) == 1 else result

        if isinstance(index, (set, frozenset)):
            index = list(index)

        result, index = [], np.atleast_1d(index)
        for table in self.loops:
            matches = index[np.any(index[:, None] == table.dtype.names, axis=1)]
            if len(matches) == 0:
                continue

            result.append(
                structured_to_unstructured(
                    table[matches], copy=True, casting="safe"
                ).squeeze(axis=1)
            )
        return _flatten_or_none(result)

    def read_cell_params(self, degrees: bool = True, normalize: bool = False):
        r"""Read the `unit cell parameters`_ (lengths and angles).

        .. _`unit cell parameters`: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ccell.html

        Parameters
        ----------
            degrees : bool, optional
                When ``True``, angles are returned in degrees (as in the CIF spec). When
                ``False``, angles are converted to radians. Default value = ``True``
            normalize: bool, optional
                Whether to scale the unit cell such that the smallest lattice parameter
                is ``1.0``.
                Default value = ``False``

        Returns
        -------
            tuple[float]:
                The box vector lengths (in angstroms) and angles (in degrees or radians)
                :math:`(L_1, L_2, L_3, \alpha, \beta, \gamma)`.

        Raises
        ------
        ValueError
            If the stored data cannot form a valid box.
        """
        angle_keys = ("_cell?angle_alpha", "_cell?angle_beta", "_cell?angle_gamma")
        box_keys = ("_cell?length_a", "_cell?length_b", "_cell?length_c", *angle_keys)

        if self.cast_values:
            cell_data = np.asarray([float(x) for x in self[box_keys]])
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
        if normalize:
            cell_data[:3] /= cell_data[:3].min()

        return tuple(float(v) for v in cell_data)  # Return as base python types

    def build_unit_cell(
        self,
        n_decimal_places: int = 4,
        additional_columns: str | Iterable[str] | None = None,
        parse_mode: str = "python_float",
        verbose: bool = False,
    ):
        """Reconstruct fractional atomic positions from Wyckoff sites and symops.

        Rather than storing an entire unit cell's atomic positions, CIF files instead
        include the data required to recreate those positions based on symmetry rules.
        Symmetry operations (stored as strings of `x,y,z` position permutations) are
        applied to the Wyckoff (symmetry irreducible) positions to create a list of
        possible atomic sites. These are then wrapped into the unit cell and filtered
        for uniqueness to yield the final crystal.

        .. tip::

            If the parsed unit cell has more atoms than expected, decrease
            ``n_decimal_places`` to account for noise. If the unit cell has fewer atoms
            than expected, increase ``n_decimal_places`` to ensure atoms are compared
            with sufficient precision. In many cases, setting ``parse_mode='sympy'``
            can improve the accuracy of reconstructed unit cells.

        Example
        -------
        Construct the atomic positions of the FCC unit cell from its Wyckoff sites:

        >>> pos = cif.build_unit_cell()
        >>> pos
        array([[0. , 0. , 0. ],
               [0. , 0.5, 0.5],
               [0.5, 0. , 0.5],
               [0.5, 0.5, 0. ]])

        Reconstruct a unit cell with its associated atomic labels. The ordering of the
        auxiliary data array will match the ordering of the atomic positions:

        >>> data = cif.build_unit_cell(additional_columns=["_atom_site_type_symbol"])
        >>> data[0] # Chemical symbol for the atoms at each lattice site
        array([['Cu'],
               ['Cu'],
               ['Cu'],
               ['Cu']], dtype='<U12')
        >>> data[1] # Lattice positions
        array([[0. , 0. , 0. ],
               [0. , 0.5, 0.5],
               [0.5, 0. , 0.5],
               [0.5, 0.5, 0. ]])
        >>> assert (pos==data[1]).all()

        Parameters
        ----------
            n_decimal_places : int, optional
                The number of decimal places to round each position to for the
                uniqueness comparison. Ideally this should be set to the number of
                decimal places included in the CIF file, but ``3`` and ``4`` work in
                most cases. Default value = ``4``
            additional_columns : str | typing.Iterable[str] | None, optional
                A column name or list of column names from the loop containing
                the Wyckoff site positions. This data is replicated alongside the atomic
                coordinates and returned in an auxiliary array.
                Default value = ``None``
            parse_mode : {'sympy', 'python_float'}, optional
                Whether to parse lattice sites symbolically (``parse_mode='sympy'``) or
                numerically (``parse_mode='python_float'``). Sympy is typically more
                accurate, but may be slower. Default value = ``'python_float'``
            verbose : bool, optional
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
        ValueError
            If the ``additional_columns`` are not properly associated with the Wyckoff
            positions.
        ImportError
            If ``parse_mode='sympy'`` and Sympy is not installed.
        """
        if parse_mode == "sympy" and not self.__class__._SYMPY_AVAILABLE:
            raise ImportError(
                "Sympy is not available! Please set parse_mode='python_float' "
                "or install sympy."
            )
        valid_modes = {"sympy", "python_float"}
        if parse_mode not in valid_modes:
            raise ValueError(f"Parse mode '{parse_mode}' not in {valid_modes}.")

        symops = self.symops

        if additional_columns is not None:
            # Find the table of Wyckoff positions and compare to keys in additional_data
            invalid_keys = next(
                set(map(str, np.atleast_1d(additional_columns))) - set(labels)
                for labels in self.loop_labels
                if set(labels) & set(self.__class__._WYCKOFF_KEYS)
            )
            if invalid_keys:
                msg = (
                    f"Requested keys {invalid_keys} are not included in the `_atom_site"
                    "_fract_[xyz]` loop and cannot be included in the unit cell."
                )
                raise ValueError(msg)

        # Read the cell params and convert to a matrix of basis vectors
        cell = self.read_cell_params(degrees=False)
        cell_matrix = _matrix_from_lengths_and_angles(*cell)

        symops_str = np.array2string(
            symops, separator=",", threshold=np.inf, floatmode="unique"
        )

        frac_strs = self.get_from_loops(self.__class__._WYCKOFF_KEYS)

        all_frac_positions = [
            _safe_eval(symops_str, *xyz, parse_mode=parse_mode) for xyz in frac_strs
        ]  # Compute N_symmetry_elements coordinates for each Wyckoff site

        pos = np.vstack(all_frac_positions)

        # Wrap into box - works generally because these are fractional coordinates
        unrounded_pos = pos.copy() % 1
        pos = pos.round(n_decimal_places) % 1

        # Filter unique points
        _, unique_fractional, unique_counts = np.unique(
            pos, return_index=True, return_counts=True, axis=0
        )

        # Double-check for duplicates with real space coordinates
        real_space_positions = pos @ cell_matrix

        _, unique_realspace, unique_counts = np.unique(
            real_space_positions.round(n_decimal_places),
            return_index=True,
            return_counts=True,
            axis=0,
        )

        # Merge unique points from realspace and fractional calculations
        unique_indices = sorted({*unique_fractional} & {*unique_realspace})

        if verbose:
            _write_debug_output(
                unique_fractional, unique_counts, pos, check="Fractional"
            )
            _write_debug_output(
                unique_fractional, unique_counts, pos, check="Realspace"
            )

        if additional_columns is None:
            return unrounded_pos[unique_indices]

        tiled_data = np.repeat(
            self.get_from_loops(additional_columns), len(symops), axis=0
        )

        return tiled_data[unique_indices], unrounded_pos[unique_indices]

    @property
    def box(self):
        """Read the unit cell as a `freud`_ or HOOMD `box-like`_ object.

        .. _`box-like`: https://hoomd-blue.readthedocs.io/en/v5.0.0/hoomd/module-box.html#hoomd.box.box_like
        .. _`freud`: https://freud.readthedocs.io/en/latest/gettingstarted/examples/module_intros/box.Box.html

        .. important::

            ``cif.box`` returns box extents and tilt factors, while
            ``CifFile.read_cell_params`` returns unit cell vector lengths and angles.
            See the `box-like`_ documentation linked above for more details.

        Example
        -------
        This method provides a convenient interface to create box objects.

        >>> box = cif.box
        >>> print(box)
        (3.6, 3.6, 3.6, 0.0, 0.0, 0.0)
        >>> import freud, hoomd # doctest: +SKIP
        >>> freud.Box(*box) # doctest: +SKIP
        freud.box.Box(Lx=3.6, Ly=3.6, Lz=3.6, xy=0, xz=0, yz=0, ...)
        >>> hoomd.Box(*box) # doctest: +SKIP
        hoomd.box.Box(Lx=3.6, Ly=3.6, Lz=3.6, xy=0.0, xz=0.0, yz=0.0)


        Returns
        -------
        tuple[float]:
            The box vector lengths (in angstroms) and unitless tilt factors.
            :math:`(L_1, L_2, L_3, xy, xz, yz)`.
        """
        return _box_from_lengths_and_angles(*self.read_cell_params(degrees=False))

    @property
    def lattice_vectors(self):
        r"""The lattice vectors of the unit cell, with :math:`\vec{a_1}\perp[100]`.

        .. important::
            The lattice vectors are stored as *columns* of the returned matrix, similar
            to `freud to_matrix()`_. This matrix must be transposed when creating a
            Freud box or transforming fractional coordinates to absolute.

        .. _`freud to_matrix()`: https://freud.readthedocs.io/en/latest/modules/box.html#freud.box.Box.to_matrix

        Example
        -------
        The box matrix can be used to transform fractional coordinates to absolute
        coordinates *after transposing to row-major form.*

        >>> lattice_vectors = cif.lattice_vectors
        >>> lattice_vectors
        array([[3.6, 0.0, 0.0],
               [0.0, 3.6, 0.0],
               [0.0, 0.0, 3.6]])
        >>> cif.build_unit_cell() @ lattice_vectors.T # Calculate absolute positions
        array([[0.0, 0.0, 0.0],
               [0.0, 1.8, 1.8],
               [1.8, 0.0, 1.8],
               [1.8, 1.8, 0.0]])

        Returns
        -------
        :math:`(3, 3)` :class:`numpy.ndarray`:
            The lattice vectors of the unit cell :math:`\vec{a_1}, \vec{a_2},\vec{a_3}`.
        """
        lx, ly, lz, xy, xz, yz = self.box
        return np.asarray([[lx, xy * ly, xz * lz], [0, ly, lz * yz], [0, 0, lz]])

    @property
    def loop_labels(self):
        """A list of column labels for each data array.

        This property is equivalent to :code:`[arr.dtype.names for arr in self.loops]`.

        Returns
        -------
        list[list[str]]:
            Column labels for :attr:`~.loops`, stored as a nested list of strings.
        """
        return [arr.dtype.names for arr in self.loops]

    @property
    def symops(self):
        r"""Extract the symmetry operations in a `parsable algebraic form`_.

        Example
        -------
        >>> cif.symops
        array([['x,y,z'],
               ['z,y+1/2,x+1/2'],
               ['z+1/2,-y,x+1/2'],
               ['z+1/2,y+1/2,x']], dtype='<U14')

        Returns
        -------
            :math:`(N,1)` numpy.ndarray[str]:
                An array containing the symmetry operations.

        .. _`parsable algebraic form`: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_symop_operation_xyz.html
        """
        # Only one key is valid in each standard, so we only ever get one match.
        return self.get_from_loops(self.__class__._SYMOP_KEYS)

    @property
    def wyckoff_positions(self):
        r"""Extract symmetry-irreducible, fractional `x,y,z` coordinates.

        Returns
        -------
            :math:`(N, 3)` :class:`numpy.ndarray`:
                Symmetry-irreducible positions of atoms in `fractional coordinates`_.

        .. _`fractional coordinates`: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html
        """
        # TODO: add additional checking in get_from_loops to verify correctness
        return cast_array_to_float(
            arr=self.get_from_loops(self.__class__._WYCKOFF_KEYS), dtype=float
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

    def _parse(self, data_iter: peekable):
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
            block = re.match(self._cpat["block_delimiter"], line.lower())
            if block is not None:
                continue

            # Extract key-value pairs and save to the internal state ===================
            pair = self._cpat["key_value_general"].match(line)

            # If we have a COD-style _key\n'long_value'
            if (
                pair is None
                and data_iter.peek("").lstrip()[:1] in "'\""
                and data_iter.peek(None)
            ):
                pair = self._cpat["key_value_general"].match(line + data_iter.peek(""))
                next(data_iter)
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
            loop = re.match(self._cpat["loop_delimiter"], line.lower())

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
                    parsed_line = [m for m in parsed_line if m != "" and m != ","]
                    loop_data.extend([parsed_line] if parsed_line else [])

                n_elements, n_cols = (
                    sum(len(row) for row in loop_data),
                    len(loop_keys),
                )

                if n_cols == 0:
                    continue  # Skip empty tables

                if n_elements % n_cols != 0:
                    warnings.warn(
                        f"Parsed data for table {len(self.loops) + 1} cannot be"
                        f" resolved into a table of the expected size and will be"
                        f"ignored. Got n={n_elements} items, expected c={n_cols}"
                        f"columns: n%c={n_elements % n_cols}).",
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
                try:
                    rectable = np.atleast_2d(loop_data)
                except ValueError as e:
                    msg = (
                        (
                            "Ragged array identified: please check the loops' syntax."
                            f"\n  Loop keys:      {loop_keys}"
                            f"\n  Processed data: {loop_data}"
                        )
                        if "setting an array element with a sequence" in str(e)
                        else e
                    )
                    raise ValueError(msg) from e
                labeled_type = [*zip(loop_keys, [dt] * n_cols)]
                try:
                    rectable.dtype = labeled_type
                except ValueError as e:
                    msg = (
                        "Loop labels do not match the structure of parsed data.\n"
                        f"  loop_labels:   {labeled_type}\n"
                        "  data[:3, ...]: "
                        f"{np.array2string(rectable[:2, :], prefix=' ' * 17)}\n"
                    )
                    raise ValueError(msg) from e
                rectable = rectable.reshape(rectable.shape, order="F")
                self.loops.append(rectable)

            if data_iter.peek(None) is None:
                break

    def __repr__(self):
        n_pairs = len(self.pairs)
        n_tabs = len(self.loops)
        return f"CifFile(fn={self._fn}) : {n_pairs} data entries, {n_tabs} data loops"

    PATTERNS: ClassVar = {
        "key_value_general": rf"^(_{_CIF_KEY}+?)\s{_PROG_PLUS}([^#]{_PROG_PLUS})",
        "loop_delimiter": rf"(loop_){_WHITESPACE}{_PROG_STAR}([^\n]{_PROG_STAR})",
        "block_delimiter": rf"(data_){_WHITESPACE}{_PROG_STAR}([^\n]{_PROG_STAR})",
        "key_list": rf"_{_CIF_KEY}{_PROG_STAR}",
        "space_delimited_data": (
            "("
            rf";[^;]*?;|"  # Non-semicolon data bracketed by semicolons
            r"'(?:'\S|[^'])*'|"  # Data with single quotes not followed by \s
            # rf"\"[^\"]{_PROG_STAR}\"|"  # Data with double quotes
            rf"[^';\"\s]{_PROG_STAR}"  # Additional non-bracketed data
            ")"
        ),
    }
    """Regex patterns used when parsing files.

    This dictionary can be modified to change parsing behavior, although doing is not
    recommended. Changes to this variable are shared across all instances of the class.
    """

    _SYMOP_KEYS = (
        "_symmetry_equiv_pos_as_xyz",
        "_space_group_symop_operation_xyz",
    )
    _WYCKOFF_KEYS = (
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_Cartn_x",
        "_atom_site_Cartn_y",
        "_atom_site_Cartn_z",
    )  # Only one set should be stored at a time
