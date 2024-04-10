.. _quickstart:

Quickstart Tutorial
===================

Once you have :ref:`installed <installation>` **parsnip**, most workflows involve reading a CIF file.
Let's assume we have the file my_file.cif in the current directory, and these are its contents:

.. literalinclude:: example_file.cif

Reading Keys
^^^^^^^^^^^^


Now, let's read extract the key-value pairs:

.. code-block:: python

    from parsnip import parse
    filename = "my_file.cif"
    pairs = parse.read_key_value_pairs(filename)
    print(pairs)
    ...    {
    ...      '_journal_year': '1999',
    ...      '_journal_page_first': '0',
    ...      '_journal_page_last': '123',
    ...      '_chemical_name_mineral': "'Copper FCC'",
    ...      '_chemical_formula_sum': "'Cu'",
    ...      '_cell_length_a': '3.6',
    ...      '_cell_length_b': '3.6',
    ...      '_cell_length_c': '3.6',
    ...      '_cell_angle_alpha': '90.0',
    ...      '_cell_angle_beta': '90.0',
    ...      '_cell_angle_gamma': '90.0'
    ...      '_symmetry_space_group_name_H-M':  'Fm-3m'
    ...    }

By default, read_key_value_pairs reads every key. To read only numeric data values, set
``only_read_numerics`` to ``True``.To take a subset, provide a tuple of strings to the ``keys`` argument.

.. code-block:: python

    # Only read the numeric data values
    pairs = parse.read_key_value_pairs(filename,only_read_numerics=True)
    print(pairs)
    ...    {
    ...      '_journal_year': 1999,
    ...      '_journal_page_first': 0,
    ...      '_journal_page_last': 123,
    ...      '_cell_length_a': 3.6,
    ...      '_cell_length_b': 3.6,
    ...      '_cell_length_c': 3.6,
    ...      '_cell_angle_alpha': 90.0,
    ...      '_cell_angle_beta': 90.0,
    ...      '_cell_angle_gamma': 90.0
    ...    }

    # Read only these keys
    keys = (
      "_journal_year"
      "_journal_page_first"
      "_journal_page_last"
    )
    pairs = parse.read_key_value_pairs(filename,keys=keys)
    print(pairs)
    ...    {
    ...      '_journal_year': '1999',
    ...      '_journal_page_first': '0',
    ...      '_journal_page_last': '123',
    ...    }

Reading Tables
^^^^^^^^^^^^^^

Now, let's read a table. To do this, we need a list of keys:

.. code-block:: python

    keys = (
      "_atom_site_label",
      "_atom_site_fract_x",
      "_atom_site_fract_y",
      "_atom_site_fract_z",
      "_atom_site_type_symbol",
      "_atom_site_Wyckoff_label"
    )
    table = parse.read_table(filename,keys=keys)
    print(table)
    ...    array([['Cu1',
    ...            '0.0000000000(0)',
    ...            '0.0000000000(0)',
    ...            '0.0000000000(0)',
    ...            'Cu'
    ...            'a']],
    ...            dtype='<U12')


Now, maybe don't need the atom site or Wyckoff labels - let's select just the numeric values, and export them as floats:

.. code-block:: python

    keys = (
      "_atom_site_fract_x",
      "_atom_site_fract_y",
      "_atom_site_fract_z",
    )
    table = parse.read_table(filename,keys=keys,cast_to_float=True)
    print(table)
    ...    array([[0., 0., 0.]], dtype=float32)

The cast_to_float argument automatically converts numeric data types, and removes tolerance and precision markers for us.
Extracting the fractional coordinates of a unit cell is a pretty common operation, so we have a convenience function that does this as well.

.. code-block:: python


    table = parse.read_fractional_positions(filename)
    print(table)
    ...    array([[0., 0., 0.]], dtype=float32)
