Initializing Molecular Simulations with ``parsnip``
===================================================

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
