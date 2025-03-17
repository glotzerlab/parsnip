Changelog
=========

The format is based on `Keep a Changelog <http://keepachangelog.com/en/1.1.0/>`__.
This project adheres to `Semantic Versioning <http://semver.org/spec/v2.0.0.html>`__.

v0.X.X - 20XX-XX-XX
-------------------

Added
~~~~~
- Additional testpath flag to conftest

Fixed
~~~~~
- Accessing data pairs with ``get_from_pairs`` or ``__getitem__`` now allows for case-insensitive searches
- Quote-delimited strings containing the delimiting character are now parsed properly

v0.2.1 - 2025-03-12
-------------------

Added
~~~~~
- New ``additional_columns`` parameter for ``build_unit_cell`` that allows the return of
  atom site labels and similar data alongside unit cell positions.
- Ensured consistent ordering of lattice positions returned from ``build_unit_cell``.
- CI testing on Windows and macOS

Fixed
~~~~~
- Type hints now properly link to their associated documentation.

v0.2.0 - 2025-02-19
-------------------

Added
~~~~~
- Support for nonsimple (';'-delimited) data entries.
- Improved support for entries containing special characters.
- Ability to query multiple keys or columns simultaneously.
- Additional tests for AMCSD and zeolite databases.
- Additional documentation and examples for the new interface

Changed
~~~~~~~
- Primary interface is now the ``CifFile`` object, which supports all previously implemented features in addition to several new methods.
- Files are now parsed lazily, and are traversed a single time.

Dependencies
~~~~~~~~~~~~
- Added ``more-itertools`` as a dependency for ``peekable`` iterators


v0.1.0 - 2024-12-20
-------------------

Added
~~~~~
- Unitcells module
- Function-based parsing interface for key and table reading
