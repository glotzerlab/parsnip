.. _quickstart:

Quickstart Tutorial
===================

Once you have :ref:`installed <installation>` **parsnip**, most workflows involve reading a CIF file.
Let's assume we have the file my_file.cif in the current directory, and these are its contents:

.. literalinclude:: example_file.cif

Reading Keys
^^^^^^^^^^^^

Now, let's read extract the key-value pairs from our cif file. This subset of data
usually contains information to reconstruct the system's unit cell, and provides
information regarding the origin of the data.


.. testsetup::

    >>> import os
    >>> if "doc/source" not in os.getcwd(): os.chdir("doc/source")

.. doctest::

    >>> from parsnip import CifFile
    >>> filename = "example_file.cif"
    >>> cif = CifFile(filename)

    >>> cif.pairs
    {'_journal_year': '1999',
    '_journal_page_first': '0',
    '_journal_page_last': '123',
    '_chemical_name_mineral': "'Copper FCC'",
    '_chemical_formula_sum': "'Cu'",
    '_cell_length_a': '3.6',
    '_cell_length_b': '3.6',
    '_cell_length_c': '3.6',
    '_cell_angle_alpha': '90.0',
    '_cell_angle_beta': '90.0',
    '_cell_angle_gamma': '90.0',
    '_symmetry_space_group_name_H-M': "'Fm-3m'"}


A `dict`-like getter syntax is provided to key-value pairs. Single keys function exactly
as a python dict, while lists of keys return lists of values. Keys not present in the
:attr:`~.pairs` dict instead return :code:`None`.

.. doctest::

    >>> cif["_journal_year"]
    '1999'

    >>> assert cif["_not_in_pairs"] is None

    # Multiple keys can be accessed simultaneously!
    >>> cif[["_cell_length_a", "_cell_length_b", "_cell_length_c"]]
    ['3.6', '3.6', '3.6']

Note that all data is stored and returned as strings by default. It is not generally
feasible to determine whether a piece of data should be processed, as conversions may
be lossy. Setting the :attr:`~.cast_values` property to :code:`True` reprocesses the
data, converting to float or int where possible. Note that once data is reprocessed,
a new CifFile object must be created to restore the original string data

.. doctest::

    >>> cif.cast_values = True # Reprocess our `pairs` dict

    >>> cif["_journal_year"]
    1999

    >>> cif[["_cell_length_a", "_cell_length_b", "_cell_length_c"]]
    [3.6, 3.6, 3.6]


Reading Tables
^^^^^^^^^^^^^^

CIF files store tables in `loop\_` delimited blocks. These structures begin with a list
of column labels (in a similar format to keys like above), followed by space-delimited
data.

This segment of the table shown above contains the table data, with 6 columns and 1 row:

.. literalinclude:: example_file.cif
    :emphasize-lines: 1
    :lines: 18-25

.. _structured arrays: https://numpy.org/doc/stable/user/basics.rec.html

Now, let's read the table. `parsnip` stores data as Numpy `structured arrays`_, which
allow for a dict-like access of data columns. The :attr:`~.loops` property returns a
list of such arrays, although the :attr:`~.get_from_loops` method is often more
convenient.


.. doctest::

    >>> len(cif.loops)
    2

    >>> cif.loops[0]
    array([[('Cu1', '0.0000000000', '0.0000000000', '0.0000000000', 'Cu', 'a')]],
    dtype=[('_atom_site_label', '<U12'),
    ('_atom_site_fract_x', '<U12'),
    ('_atom_site_fract_y', '<U12'),
    ('_atom_site_fract_z', '<U12'),
    ('_atom_site_type_symbol', '<U12'),
    ('_atom_site_Wyckoff_label', '<U12')])

    >>> cif.loops[0]["_atom_site_label"]
        array([['Cu1']], dtype='<U12')

    # (Unstructured) slices of tables can be easily accessed!
    >>> xyz = cif.get_from_loops(["_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"])

    >>> xyz
        array([['0.0000000000', '0.0000000000', '0.0000000000']], dtype='<U12')

    >>> xyz.astype(float)
        array([[0., 0., 0.]])




Building Unit Cells
^^^^^^^^^^^^^^^^^^^

CIF files are commonly used to reconstruct atomic positions for a particular crystal.
While the example file shown throughout this tutorial corresponds to FCC copper, it only
contains a single atomic position, in constrast to the 4 expected for FCC's
primitive cell. `parsnip` can reconstruct tilable unit cells from symmetry operations
and symmetry-irreducible (Wyckoff) positions contained in the file.

.. literalinclude:: example_file.cif
    :lines: 25
.. literalinclude:: example_file.cif
    :lines: 29-37

Only one line is required to build a tilable unit cell! The positions returned here
are in fractional coordinates, and can be imported into tools like `freud`_ to rapidly
build out supercells. For absolute coordinates (based on cell parameters stored in the
file), set :code:`fractional=False`.

.. _`freud`: https://freud.readthedocs.io/en/latest/modules/data.html#freud.data.UnitCell

.. doctest::

    >>> pos = cif.build_unit_cell(fractional=True)
    >>> print(pos)
    [[0.  0.  0. ]
     [0.  0.5 0.5]
     [0.5 0.  0.5]
     [0.5 0.5 0. ]]

Once `freud`_ is installed, crystal structures can be easily replicated!

.. doctest-requires:: freud

    >>> import freud
    >>> box = freud.Box(*cif.cell)
    >>> uc = freud.data.UnitCell(box, basis_positions=pos)
    >>> box, pos = uc.generate_system(num_replicas=2)

    >>> assert len(pos) == 4 * 2**3
    >>> np.testing.assert_allclose(box.L / 2, cif.cell[:3], atol=1e-15)
