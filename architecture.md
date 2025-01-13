Parsnip Architecture
--------------------

The primary design goal of ``parsnip`` was to create a lightweight, simple CIF parsing
library. ``parsnip`` has a minimal set of dependencies, relatively few lines of code,
and extensive testing to validate the accuracy of read files. Dozens of CIF parsing
libraries exist, but most are either (1) part of a much larger project (and therefore
undesirable as a simple dependency) or (2) have a poorly documented interface. This
project is designed to bridge that gap with careful documentation and a minimal subset
of features that build well into other open-source projects.


This project takes a (reasonably) permissive view of the CIF specification: data entries
that "look like" valid data will be parsed, regardless of file encoding, line length,
special characters, or syntax specifics like unlabeled blocks. ``parsnip`` is designed
to read and extract data from CIF, mmCIF, and STAR files regardless of their compliance
with the full CIF specification.
