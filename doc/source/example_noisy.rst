Reconstrucing Noisy Unit Cells
==============================

Diffraction experiments and other experimental techniques for quantifying structure
typically offer limited precision in the measurements that can be made. As a result,
the Wyckoff position data recorded in some CIF files -- particularly older ones -- may
make reproduction of the original structure challenging. In this example, we explore how
**parsnip**'s `build_unit_cell` method can be tuned to accurately reproduce structures
with complicated geometries, using alpha-Selenium as an example.

.. testsetup::

    >>> import os
    >>> import numpy as np
    >>> if "doc/source" not in os.getcwd(): os.chdir("doc/source")


.. literalinclude:: hP3.cif

Note that the basis positions for alpha-Selenium are provided to five decimal
places of accuracy, while the symmetry operations are provided in a rational form.

.. doctest::

    >>> from parsnip import CifFile
    >>> cif = CifFile("hP3.cif")
    >>> # Let's make sure we reconstruct the unit cell's three atoms
    >>> correct_uc = cif.build_unit_cell()
    >>> correct_uc
    array([[0.2254    , 0.        , 0.33333   ],
           [0.        , 0.2254    , 0.66666333],
           [0.7746    , 0.7746    , 0.99999667]])
    >>> site_multiplicity = int(cif["_atom_site_symmetry_multiplicity"].squeeze())
    >>> assert len(correct_uc) == site_multiplicity

**parsnip**'s default settings are able to correctly reproduce the unit cell -- but
the mismatch between numerical data and the symmetry operation strings can cause issues.
If we truncate the Wyckoff position data, even by one decimal place, the reconstructed
crystal contains duplicate atoms:

.. literalinclude:: hP3-four-decimal-places.cif
   :diff: hP3.cif

Rebuilding our crystal results in an error:

.. doctest::

    >>> lower_precision_cif = CifFile("hP3-four-decimal-places.cif")
    >>> uc = lower_precision_cif.build_unit_cell(n_decimal_places=4)
    >>> uc # doctest: +SKIP
    array([[0.2254    , 0.        , 0.3333    ],  # A
           [0.        , 0.2254    , 0.66663333],  # B
           [0.7746    , 0.7746    , 0.99996667],  # C
           [0.2254    , 0.        , 0.33336667],  # A
           [0.        , 0.2254    , 0.6667    ]]) # B
    >>> uc.shape == correct_uc.shape # Our unit cell has duplicate atoms!
    False

By default, **parsnip** uses three decimal places of accuracy to reconstruct crystals.
This yields the best overall accuracy (tested with several thousand CIFs), but is not
always the best choice in general. A good rule of thumb is to use one fewer decimal
places than the CIF file contains. This ensures positions are rounded sufficiently to
detect duplicate atoms, and avoids issues in complex structures where Wyckoff positions
may be very close to one another. Using this setting results in the correct structure
once again.


.. doctest::

    >>> cif = CifFile("hP3-four-decimal-places.cif")
    >>> three_decimal_places = cif.build_unit_cell(n_decimal_places=3)
    >>> three_decimal_places
    array([[0.2254    , 0.        , 0.3333    ],
           [0.        , 0.2254    , 0.66663333],
           [0.7746    , 0.7746    , 0.99996667]])
    >>> assert three_decimal_places.shape == correct_uc.shape

.. important::

    Rounding of Wyckoff positions is an intermediate step in the unit cell
    reconstruction, and does not negatively impact the accuracy of the returned data.
    The unit cell is always returned in the full precision of the input CIF:

    .. doctest::

        >>> cif = CifFile("hP3-four-decimal-places.cif")
        >>> one_decimal_place = cif.build_unit_cell(n_decimal_places=1)
        >>> np.testing.assert_array_equal(one_decimal_place, three_decimal_places)


Increasing Performance
^^^^^^^^^^^^^^^^^^^^^^

In some cases, particularly when constructing thousands of unti cells, the performance
of parsnip's ``build_unit_cell`` may become a bottleneck. **parsnip** includes several
tools for resolving this: first, ``parse_mode="python_float"`` attempts to build unit
cells using floating point arithmetic rather than rational expression. This is less
accurate, but is still sufficient for high-quality databases and stuctures. For the
best combination of performance and accuracy, installing the `cfractions`_ library lets
**parsnip** use more optimized code for unit cell reconstruction.

.. _cfractions: https://pypi.org/project/cfractions/

.. doctest::

    >>> # uv pip install cfractions
    >>> cif = CifFile("hP3.cif")
    >>> faster = cif.build_unit_cell(n_decimal_places=4)
    >>> faster
    array([[0.2254    , 0.        , 0.33333  ],
           [0.        , 0.2254    , 0.6666633],
           [0.7746    , 0.7746    , 0.99999667]])
    >>> assert faster.shape == correct_uc.shape

.. _sympy: https://www.sympy.org/en/index.html
