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
    array([[0.2254    , 0.        , 0.33333333],
           [0.        , 0.2254    , 0.66666667],
           [0.7746    , 0.7746    , 0.        ]])
    >>> site_multiplicity = int(cif["_atom_site_symmetry_multiplicity"].squeeze())
    >>> assert len(correct_uc) == site_multiplicity

**parsnip**'s default settings are able to correctly reproduce the unit cell -- but
the mismatch between numerical data and the symmetry operation strings can cause issues.
If we truncate the Wyckoff position data, even by one decimal place, we might expect
the reconstructed crystal to contain duplicate atoms:

.. literalinclude:: hP3-four-decimal-places.cif
   :diff: hP3.cif

However, **parsnip**'s rational mode recognizes truncated representations of exact
fractions (e.g., ``0.3333`` is identified as ``1/3``) and snaps them to their exact
values before applying symmetry operations. This means the truncated file still
produces the correct structure:

.. doctest::

    >>> lower_precision_cif = CifFile("hP3-four-decimal-places.cif")
    >>> uc = lower_precision_cif.build_unit_cell(n_decimal_places=4)
    >>> uc
    array([[0.2254    , 0.        , 0.33333333],
           [0.        , 0.2254    , 0.66666667],
           [0.7746    , 0.7746    , 0.        ]])
    >>> uc.shape == correct_uc.shape
    True

For structures where coordinates are not near ideal fractions, the ``n_decimal_places``
parameter controls deduplication precision. By default, **parsnip** uses three decimal
places, which yields the best overall accuracy. A good rule of thumb is to use one
fewer decimal place than the CIF file contains, ensuring positions are rounded
sufficiently to detect duplicate atoms in complex structures.

.. important::

    Rounding of Wyckoff positions is an intermediate step in the unit cell
    reconstruction, and does not negatively impact the accuracy of the returned data.
    The unit cell is always returned in the full precision of the input CIF:

    .. doctest::

        >>> cif = CifFile("hP3-four-decimal-places.cif")
        >>> one_decimal_place = cif.build_unit_cell(n_decimal_places=1)
        >>> np.testing.assert_array_equal(one_decimal_place, uc)


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
    array([[0.2254    , 0.        , 0.33333333],
           [0.        , 0.2254    , 0.66666667],
           [0.7746    , 0.7746    , 0.        ]])
    >>> assert faster.shape == correct_uc.shape

.. _sympy: https://www.sympy.org/en/index.html
