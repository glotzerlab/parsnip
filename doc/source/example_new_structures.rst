.. _setbasis:

Refining and Experimenting with Structures
==========================================

**parsnip** allows users to set the Wyckoff positions of a crystal, enabling the
construction of modified (or entirely new) structures. In this example, we show
how an experimental β-Manganese structure can be refined into the more uniform variant
described by `O'Keefe and Andersson`_.

.. _`O'Keefe and Andersson`: https://doi.org/10.1107/S0567739477002228


These are the Wyckoff positions for elemental β-Manganese, drawn directly from a CIF
file:

.. literalinclude:: betamn.cif
   :lines: 51-52

Formatted more nicely, we see the following:

.. list-table:: Crystallographic data loop for β-Mn
   :widths: 15 15 20 15 10 10 10
   :header-rows: 1

   * - Site Label
     - Type Symbol
     - Symmetry Multiplicity
     - Wyckoff letter
     - x
     - y
     - z
   * - Mn1
     - Mn
     - 8
     - c
     - 0.06361
     - 0.06361
     - 0.06361
   * - Mn2
     - Mn
     - 12
     - d
     - 0.12500
     - 0.20224
     - 0.45224

First, we note the symmetry multiplicity (8 for :math:`Mn_1` and 12 for
:math:`Mn_2`), which indicates how many atomic positions arise from each Wyckoff
site. Second, we can identify that each Wyckoff position is labeled by a letter that
differentiates it from other sites. While this tutorial will not delve too deeply into
crystallography, it is sufficient to note that this Wyckoff letter provides a mapping
to the `International Tables`_ for each space group. For β-Manganese, we will use this
mapping to identify one coordinate equation that describes each site. The correct table
for β-Mn is included in the tabs below, with the coordinate equations that match the
CIF data (:math:`(x, x, x)` and :math:`(1/8, y, y + 1/4)`) highlighted in bold on their
corresponding tabs.

.. _`International Tables`: https://web.archive.org/web/20170430110556/http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list?gnum=213

.. tab:: 4a

   .. list-table::
      :widths: 30 70
      :header-rows: 0

      * - **Multiplicity**
        - 4
      * - **Site Symmetry**
        - ``.32``
      * - **Coordinates**
        - .. list-table::
             :header-rows: 0
             :widths: 50 50

             * - :math:`(3/8, 3/8, 3/8)`
               - :math:`(1/8, 5/8, 7/8)`
             * - :math:`(5/8, 7/8, 1/8)`
               - :math:`(7/8, 1/8, 5/8)`

.. tab:: 4b

   .. list-table::
      :widths: 30 70
      :header-rows: 0

      * - **Multiplicity**
        - 4
      * - **Site Symmetry**
        - ``.32``
      * - **Coordinates**
        - .. list-table::
             :header-rows: 0
             :widths: 50 50

             * - :math:`(7/8, 7/8, 7/8)`
               - :math:`(5/8, 1/8, 3/8)`
             * - :math:`(1/8, 3/8, 5/8)`
               - :math:`(3/8, 5/8, 1/8)`

.. tab:: 8c

   .. list-table::
      :widths: 30 70
      :header-rows: 0

      * - **Multiplicity**
        - 8
      * - **Site Symmetry**
        - ``.3.``
      * - **Coordinates**
        - .. list-table::
             :header-rows: 0
             :widths: 50 50

             * - :math:`\mathbf{(x, x, x)}`
               - :math:`(-x+1/2, -x, x+1/2)`
             * - :math:`(-x, x+1/2, -x+1/2)`
               - :math:`(x+1/2, -x+1/2, -x)`
             * - :math:`(x+3/4, x+1/4, -x+1/4)`
               - :math:`(-x+3/4, -x+3/4, -x+3/4)`
             * - :math:`(x+1/4, -x+1/4, x+3/4)`
               - :math:`(-x+1/4, x+3/4, x+1/4)`

.. tab:: 12d

   .. list-table::
      :widths: 30 70
      :header-rows: 0

      * - **Multiplicity**
        - 12
      * - **Site Symmetry**
        - ``..2``
      * - **Coordinates**
        - .. list-table::
             :header-rows: 0
             :widths: 50 50

             * - :math:`\mathbf{(1/8, y, y+1/4)}`
               - :math:`(3/8, -y, y+3/4)`
             * - :math:`(7/8, y+1/2, -y+1/4)`
               - :math:`(5/8, -y+1/2, -y+3/4)`
             * - :math:`(y+1/4, 1/8, y)`
               - :math:`(y+3/4, 3/8, -y)`
             * - :math:`(-y+1/4, 7/8, y+1/2)`
               - :math:`(-y+3/4, 5/8, -y+1/2)`
             * - :math:`(y, y+1/4, 1/8)`
               - :math:`(-y, y+3/4, 3/8)`
             * - :math:`(y+1/2, -y+1/4, 7/8)`
               - :math:`(-y+1/2, -y+3/4, 5/8)`

.. tab:: 24e

   .. list-table::
      :widths: 30 70
      :header-rows: 0

      * - **Multiplicity**
        - 24
      * - **Site Symmetry**
        - ``1``
      * - **Coordinates**
        - .. list-table::
             :header-rows: 0
             :widths: 50 50

             * - :math:`(x, y, z)`
               - :math:`(-x+1/2, -y, z+1/2)`
             * - :math:`(-x, y+1/2, -z+1/2)`
               - :math:`(x+1/2, -y+1/2, -z)`
             * - :math:`(z, x, y)`
               - :math:`(z+1/2, -x+1/2, -y)`
             * - :math:`(-z+1/2, -x, y+1/2)`
               - :math:`(-z, x+1/2, -y+1/2)`
             * - :math:`(y, z, x)`
               - :math:`(-y, z+1/2, -x+1/2)`
             * - :math:`(y+1/2, -z+1/2, -x)`
               - :math:`(-y+1/2, -z, x+1/2)`
             * - :math:`(y+3/4, x+1/4, -z+1/4)`
               - :math:`(-y+3/4, -x+3/4, -z+3/4)`
             * - :math:`(y+1/4, -x+1/4, z+3/4)`
               - :math:`(-y+1/4, x+3/4, z+1/4)`
             * - :math:`(x+3/4, z+1/4, -y+1/4)`
               - :math:`(-x+1/4, z+3/4, y+1/4)`
             * - :math:`(-x+3/4, -z+3/4, -y+3/4)`
               - :math:`(x+1/4, -z+1/4, y+3/4)`
             * - :math:`(z+3/4, y+1/4, -x+1/4)`
               - :math:`(z+1/4, -y+1/4, x+3/4)`
             * - :math:`(-z+1/4, y+3/4, x+1/4)`
               - :math:`(-z+3/4, -y+3/4, -x+3/4)`

.. testsetup::

    >>> import os
    >>> import numpy as np
    >>> if "doc/source" not in os.getcwd(): os.chdir("doc/source")

Loading the file shows the twenty atoms we expect for β-Mn:

.. doctest::

    >>> from parsnip import CifFile
    >>> filename = "betamn.cif"
    >>> cif = CifFile(filename)
    >>> uc = cif.build_unit_cell()
    >>> assert uc.shape == (8 + 12, 3)

And of course, the Wyckoff position data reflects the data tabulated above:

    >>> mn1, mn2 = cif.wyckoff_positions
    >>> mn1
    array([0.06361, 0.06361, 0.06361])
    >>> mn2
    array([0.125  , 0.20224, 0.45224])
    >>> x = mn1[0]
    >>> y = mn2[1]
    >>> np.testing.assert_allclose(mn1, x)
    >>> np.testing.assert_allclose(mn2[2], y + 1 / 4)
    >>> np.testing.assert_allclose(mn2[0], 1 / 8)

Exploring β-Manganese
^^^^^^^^^^^^^^^^^^^^^

β-Manganese is a `tetrahedrally close-packed`_ (TCP) structure, a class of complex
phases whose geometry minimizes the distance between atoms in a manner that prevents the
formation of octahedral interstitial sites. Intuitively, one can image the bond network
of TCP structures forming a space-filling collection of irregular tetrahedra, with some
required amount of distortion imposed by the requirement that the structure tiles space.

It turns out that natural β-Manganese actually has *more* variation in bond lengths
than is strictly required for this topology of structure. `O'Keefe and Andersson`_
noticed that moving the ``Mn1`` and ``Mn2`` Wyckoff positions by just ``0.0011`` and
``0.0042`` fractional units results in a TCP structure composed of only five unique bond
lengths, rather than the seven observed in experimental β-Manganese.

.. _`tetrahedrally close-packed`: https://www.chemie-biologie.uni-siegen.de/ac/hjd/lehre/ss08/vortraege/mehboob_tetrahedrally_close_packing_corr_.pdf

Using **parsnip**, we can explore the differences between experimental and ideal
β-Manganese, quantifying the distribution of bond lengths in the crystal:

.. doctest::

    >>> from parsnip import CifFile
    >>> from math import sqrt
    >>> filename = "betamn.cif"
    >>> cif = CifFile(filename)
    >>> atomic_uc = cif.build_unit_cell()
    >>> assert atomic_uc.shape == (20, 3)
    >>> # Values are drawn from O'Keefe and Andersson, linked above.
    >>> x = 1 / (9 + sqrt(33))      # Parameter for the 8c Wyckoff position
    >>> mn1 = [x, x, x] # doctest: +FLOAT_CMP
    >>> mn1
    [0.0678216, 0.0678216, 0.0678216]
    >>> y = (9 - sqrt(33)) / 16
    >>> mn2 = [1 / 8, y, y + 1 / 4] # Parameter for the 12d Wyckoff position
    >>> mn2 # doctest: +FLOAT_CMP
    [0.1250000, 0.2034648, 0.4534648]

    >>> _ = cif.set_wyckoff_positions([mn1, mn2])
    >>> # We should still have the same number of atoms
    >>> ideal_uc = cif.build_unit_cell(n_decimal_places=4)
    >>> assert ideal_uc.shape == atomic_uc.shape


Analyzing our New Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following plot shows a histogram of neighbor distances for experimental
β-Manganese (top) and the ideal structure (bottom). Each bar corresponds with a
single neighbor bond length, with each particle's neighbors existing at one of the
specified distances. Interestingly, althought the ideal structure has a more uniform
topology with fewer total distinct edges, the observed atomic structure more uniformly
distributes bonds to each particle.


.. image:: _static/perfect_imperfect_bmn.svg
  :width: 100%


A Note on Symmetry
^^^^^^^^^^^^^^^^^^

Modifying the Wyckoff positions of a crystal (without changing the symmetry operations)
cannot reduce the symmetry of the structure --- however, some choices of sites can
result in *additional* symmetry operations that are not present in the input space
group. While the example provided above preserved the space group of our crystal,
choosing a fractional coordinate that lies on a high symmetry point (like the origin,
or the center of the cell) can result in differences.


.. doctest-requires:: spglib

    >>> import spglib
    >>> box = cif.lattice_vectors
    >>> # Verify that our initial and "ideal" β-Manganese cells share a space group
    >>> spglib.get_spacegroup((box, atomic_uc, [0] * 20))
    'P4_132 (213)'
    >>> spglib.get_spacegroup((box, ideal_uc, [0] * 20))
    'P4_132 (213)'
    >>> cif["_symmetry_Int_Tables_number"] # Data from the initial file.
    '213'


Placing a Wyckoff position on a high-symmetry site results in a change in the space
group.

.. doctest-requires:: spglib

    >>> cif = CifFile("betamn.cif").set_wyckoff_positions([[0.0, 0.0, 0.0]])
    >>> different_uc = cif.build_unit_cell()
    >>> spglib.get_spacegroup((box, different_uc, [0] * len(different_uc)))
    'Fd-3m (227)'

Design Rules for Crystal Construction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While the global space group symmetry can only be increased by changing the Wyckoff
positions, the point group symmetry of sites can vary greatly. The example above
chose points that maintained the multiplicity of each site, but general choices do not
preserve this. First, let's confirm that the Wyckoff letters and site point groups are
the same in the atomic and ideal crystals:

.. doctest-requires:: spglib

    >>> def get_particle_point_groups(box, basis):
    ...     spglib_cell = (box, basis, [0] * len(basis))
    ...     dataset = spglib.get_symmetry_dataset(spglib_cell)
    ...     wycks = sorted({*dataset.wyckoffs})
    ...     point_groups = sorted({*dataset.site_symmetry_symbols})[::-1]
    ...     return (wycks, point_groups)
    >>> get_particle_point_groups(box, atomic_uc)
    (['c', 'd'], ['.3.', '..2'])
    >>> get_particle_point_groups(box, ideal_uc)
    (['c', 'd'], ['.3.', '..2'])


A more general choice of the basis will often result in different point symmetry.
Referring to the `International Tables`_ for space group 213 shows the ``a`` and ``b``
Wyckoff positions, which have higher symmetry and a lower multiplicity. Selecting any
value from the "coordinates" table for the 4a position yields the expected 4-particle
unit cell with a site symmetry of ``'.32'``.


.. doctest-requires:: spglib

    >>> four_a = [[3/8, 3/8, 3/8]]
    >>> four_a_cif = CifFile("betamn.cif").set_wyckoff_positions(four_a)
    >>> four_a_uc = four_a_cif.build_unit_cell()
    >>> get_particle_point_groups(box, four_a_uc)
    (['a'], ['.32'])
    >>> assert four_a_uc.shape == (4, 3)
    >>> spglib.get_spacegroup((box, four_a_uc, [0] * 4))
    'P4_132 (213)'
    >>> four_a_uc
    array([[0.375, 0.375, 0.375],
           [0.875, 0.125, 0.625],
           [0.625, 0.875, 0.125],
           [0.125, 0.625, 0.875]])


When working with systems where the same particle type lies on multiple Wyckoff
positions, care must be taken to ensure those sites do not satisfy a symmetry operation
in a *higher* space group than the target. The following example assigns the 4a and 4b
Wyckoff positions to a single atomic type. Even though the reconstructed crystal
contains the expected 8 particles, the sites are related by the symmetry element
``x+1/2, y+1/2, z+1/2`` of the next highest space group, #214.


.. doctest-requires:: spglib

    >>> four_a_four_b = [[3/8, 3/8, 3/8], [7/8, 7/8, 7/8]]
    >>> four_a_four_b_cif = CifFile("betamn.cif").set_wyckoff_positions(four_a_four_b)
    >>> four_a_four_b_uc = four_a_four_b_cif.build_unit_cell()
    >>> assert four_a_four_b_uc.shape == (8, 3)
    >>> # NOTE: these sites are equivalent under a *higher* space group!
    >>> get_particle_point_groups(box, four_a_four_b_uc)
    (['b'], ['.32'])
    >>> spglib.get_spacegroup((box, four_a_four_b_uc, [0] * 8))
    'I4_132 (214)'
    >>> # If the sites are different elements, the space group is preserved
    >>> spglib.get_spacegroup((box, four_a_four_b_uc, [0,0,0,0, 1,1,1,1]))
    'P4_132 (213)'

A similar consideration must be made for Wyckoff positions whose coordinates contain
one or more degrees of freedom. In β-Manganese, the 8c and 12d Wyckoff sites each
have one degree of freedom -- the ``x`` and ``y`` variables assigned above. If we set
these degrees of freedom such that Wyckoff positions are no longer independent, we also
alter the space group of the structure. In this case, we solve the system of equations
that arises from setting the coordinates :math:`[x, x, x] = [1 / 8, y, y + 1 / 4]` and
assign that value to both ``x`` and ``y``. The resulting points end up reconstructing
the 16d Wyckoff position in the space group #227!

.. doctest-requires:: spglib

    >>> x = y = -1/8
    >>> wyckoff_c = [x, x, x]
    >>> wyckoff_d = [1 / 8, y, y + 1 / 4]
    >>> c_d_linked = [wyckoff_c, wyckoff_d]
    >>> not_beta_manganese = CifFile("betamn.cif").set_wyckoff_positions(c_d_linked)
    >>> not_beta_mn_uc = not_beta_manganese.build_unit_cell()
    >>> not_beta_mn_uc.shape # NOTE: this is no longer 12+8 sites!
    (16, 3)
    >>> spglib.get_spacegroup((box, not_beta_mn_uc, [0] * 16))
    'Fd-3m (227)'
    >>> get_particle_point_groups(box, not_beta_mn_uc)
    (['d'], ['.-3m'])

.. tab:: #227 16d

   .. list-table::
      :widths: 30 70
      :header-rows: 0

      * - **Space Group**
        - #227
      * - **Multiplicity**
        - 16
      * - **Site Symmetry**
        - ``.-3m``
      * - **Coordinates**
        - .. list-table::
             :header-rows: 0
             :widths: 50 50

             * - :math:`(5/8, 5/8, 5/8)`
               - :math:`(3/8, 7/8, 1/8)`
             * - :math:`(7/8, 1/8, 3/8)`
               - :math:`(1/8, 3/8, 7/8)`

Takeaways
^^^^^^^^^

The examples above give rise to a few design rules for structure refinement and
modification:

1. For Wyckoff positions without degrees of freedom, care must be taken to ensure sites are not linked by symmetry operations present in a higher space group.
    - In general, this can be verified be comparing against a list of operations from
      the IUCR Crystal database, or a CIF file with space group #230.
2. For Wyckoff positions *with* degrees of freedom, the following must be ensured:
    - Wyckoff positions with different labels must be linearly independent (i.e. their
      coordinate equations must not be equal for the chosen degrees of freedom).
    - Free variables must be chosen such that the points do not lie on high-symmetry
      locations, particularly the origin and power-of-two fractions. This condition is
      equivalent to that in point (1), and may be resolved in a similar manner.


**parsnip** allows us to use existing structural data to generate new crystals,
including those that have not been observed in experiment. While the example shown here
is relatively simple, assigning alternative Wyckoff positions enables high-throughput
materials discovery research and offers a simple framework by which structural features
can be explored.
