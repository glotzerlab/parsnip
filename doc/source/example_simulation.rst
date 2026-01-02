Initializing Molecular Simulations
==================================

When performing molecular simulations of solid materials, it is often useful to
initialize your system in a crystal structure. **parsnip** makes this extremely easy!

HOOMD-Blue
^^^^^^^^^^

HOOMD-Blue can operate directly on array data, so we can move data directly from
**parsnip** to the simulation itself.

.. testsetup::

    >>> import os
    >>> import numpy as np
    >>> if "doc/source" not in os.getcwd(): os.chdir("doc/source")
    >>> from parsnip import CifFile
    >>> filename = "example_file.cif"
    >>> cif = CifFile(filename)
    >>> # Mock HOOMD import, as it is not available via pip install
    >>> import gsd.hoomd as hoomd
    >>> hoomd.Snapshot = hoomd.Frame
    >>> snapshot = hoomd.Snapshot()
    >>> # Pre-initialize the data arrays, as gsd does not support HOOMD's arr[:]
    >>> #    pattern for assignment. This data will be overwritten in the doctest.
    >>> snapshot.particles.position = np.full((4, 3), -999.0)
    >>> snapshot.particles.typeid = np.full((4,), -999)

.. doctest::

    >>> import hoomd # doctest: +SKIP
    >>> from parsnip import CifFile
    >>> filename = "example_file.cif"
    >>> cif = CifFile(filename)

    >>> snapshot = hoomd.Snapshot() # doctest: +SKIP
    >>> snapshot.particles.N = len(cif.build_unit_cell())
    >>> snapshot.particles.position[:] = cif.build_unit_cell()
    >>> snapshot.configuration.box = cif.box
    >>> snapshot.particles.types = ["Particle"]

    >>> snapshot.replicate(nx=2, ny=2, nz=3) # 2 x 2 x 3 supercell # doctest: +SKIP
    >>> assert snapshot.particles.N == (2 * 2 * 3) * len(pos)      # doctest: +SKIP


Once the snapshot is constructed, it can be attached to a simulation as follows:

.. doctest-skip::

    >>> import hoomd
    >>> simulation = hoomd.Simulation(device=hoomd.device.CPU())
    >>> simulation.create_state_from_snapshot(snapshot)

If we want to extract additional data for our simulation, there are a few extra steps.
In HOOMD-Blue, ``particle.types`` are unique string identifiers that get mapped to
individual particles via the ``particles.typeid`` array. The following code extracts
``_atom_site_type_label`` data and assigns the "Cu" atom type to all particles. For
structures with multiple atom sites, the ``particles.typeid`` array will have nonzero
entries that correspond with other type labels.


.. doctest::

    >>> from collections import defaultdict

    >>> labels, pos = cif.build_unit_cell(additional_columns=["_atom_site_type_symbol"])
    >>> # ... initialize the snapshot's `N`, `box`, and `position` data as above

    >>> particle_type_map = defaultdict(list)
    >>> for i, label in enumerate(labels.squeeze(axis=1)):
    ...     particle_type_map[label].append(i)
    >>> particle_type_map["Cu"] # Atoms 1-4 have the type symbol "Cu"
    [0, 1, 2, 3]

    >>> # Construct the TypeIDs that map our atomic symbol to the corresponding position
    >>> typeid_array = np.ones(len(snapshot.particles.position), dtype=int)
    >>> for typeid, label in enumerate(particle_type_map.keys()):
    ...     typeid_array[particle_type_map[label]] = typeid
    >>> snapshot.particles.typeid[:] = typeid_array
    >>> snapshot.particles.typeid
    array([0, 0, 0, 0])

    >>> snapshot.particles.types = [str(key) for key in particle_type_map.keys()]
    >>> snapshot.particles.types
    ['Cu']

    >>> assert len(snapshot.particles.types) == len(cif["_atom_site_type_symbol"])


LAMMPS
^^^^^^

In contrast to HOOMD-Blue, LAMMPS typically requires us to write out structure data to
a `LAMMPS Data File`_ before simulations can begin. Although topology data is not
commonly stored in CIF files, **parsnip** makes it simple to reconstruct atomic crystals
in LAMMPS.


.. _`LAMMPS Data File`: https://docs.lammps.org/2001/data_format.html

.. doctest::

    >>> from collections import defaultdict

    >>> def write_lammps_data(cif: CifFile, atom_type_labels: bool = True):
    ...     """Convert a CIF file into a LAMMPS data file."""
    ...     data = "(LAMMPS Data File, written with parsnip)\n\n"
    ...
    ...     fractional_coordinates = cif.build_unit_cell()
    ...     atomic_positions = fractional_coordinates @ cif.lattice_vectors.T
    ...
    ...     atom_types = cif["_atom_site_type_symbol"]
    ...     particle_type_map = defaultdict(list)
    ...
    ...     # Write out the number of atoms and atom types
    ...     data += f"{len(atomic_positions)} atoms\n"
    ...     data += f"{len(atom_types)} atom types\n\n"
    ...
    ...     # Write out the box, including the (zero, in this case) tilt factors
    ...     lx, ly, lz, xy, xz, yz = cif.box
    ...     data += f"0.0 {lx:.12f} xlo xhi\n"
    ...     data += f"0.0 {ly:.12f} ylo yhi\n"
    ...     data += f"0.0 {lz:.12f} zlo zhi\n"
    ...     data += f"{xy:.12f} {xz:.12f} {yz:.12f} xy xz yz\n\n"
    ...
    ...     for i, label in enumerate(labels.squeeze(axis=1)):
    ...         particle_type_map[label].append(i)
    ...
    ...     # Write out atom type labels
    ...     if atom_type_labels:
    ...         data += "\nAtom Type Labels\n\n"
    ...         for i, unique_labels in enumerate(particle_type_map):
    ...             data += f"{i} {label}\n"
    ...         data += "\n"
    ...
    ...     # Write out the atomic position data -- note the similarities with typeid!
    ...     data += f"Atoms # atomic\n\n"
    ...
    ...     # Construct the TypeIDs that map our atomic symbol to an index
    ...     atom_typeid_array = np.ones(len(atomic_positions), dtype=int)
    ...     atom_type_array = np.ones(len(atomic_positions), dtype="U4")
    ...     for typeid, label in enumerate(particle_type_map.keys()):
    ...         atom_typeid_array[particle_type_map[label]] = typeid
    ...         atom_type_array[particle_type_map[label]] = label
    ...
    ...     for i, coordinate in enumerate(atomic_positions):
    ...         coord_str = " ".join([f"{xyz:.12f}" for xyz in coordinate])
    ...         label = atom_type_array[i] if atom_type_labels else atom_typeid_array[i]
    ...         data += f"  {i}   {label}  {coord_str}\n"
    ...
    ...     return data

    >>> print(write_lammps_data(cif))
    (LAMMPS Data File, written with parsnip)
    <BLANKLINE>
    4 atoms
    1 atom types
    <BLANKLINE>
    0.0 3.600000000000 xlo xhi
    0.0 3.600000000000 ylo yhi
    0.0 3.600000000000 zlo zhi
    0.000000000000 0.000000000000 0.000000000000 xy xz yz
    <BLANKLINE>
    Atom Type Labels
    <BLANKLINE>
    0 Cu
    <BLANKLINE>
    Atoms # atomic
      0   Cu  0.000000000000 0.000000000000 0.000000000000
      1   Cu  0.000000000000 1.800000000000 1.800000000000
      2   Cu  1.800000000000 0.000000000000 1.800000000000
      3   Cu  1.800000000000 1.800000000000 0.000000000000

.. Validate our output data is (1) valid LAMMPS data and (2) reconstructs our system.
.. testcleanup::

    >>> from io import StringIO
    >>> from packaging import version
    >>> import ase
    >>> from ase.io import read
    >>> write_type = version.parse(ase.__version__) >= version.parse("3.27.0")
    >>> atoms = read(StringIO(write_lammps_data(cif, write_type)), format='lammps-data')
    >>> fractional_coordinates = cif.build_unit_cell()
    >>> atomic_positions = fractional_coordinates @ cif.lattice_vectors.T
    >>> assert len(atomic_positions) == 4
    >>> num = 29 if write_type else 0
    >>> np.testing.assert_array_equal(atoms.get_atomic_numbers(), [num, num, num, num])
    >>> np.testing.assert_array_equal(np.diag([3.6, 3.6, 3.6]), atoms.get_cell())
