---
title: "Parsnip: Streamlined Crystallographic Data Parsing for Simulation Science"
tags:
  - Python
  - crystallography
  - materials science
authors:
  - name: Jenna Bradley
    orcid: 0009-0007-2443-2982
    affiliation: 1
  - name: Sharon Glotzer
    orcid: 0000-0002-7197-0085
    affiliation: [1, 2]
affiliations:
 - name: Materials Science and Engineering, University of Michigan, United States
   index: 1
   ror: 00jmfr291
 - name: Department of Chemical Engineering, University of Michigan
   index: 2
date: XX December 2025
bibliography: paper.bib
---

# Summary

`parsnip` provides a lightweight, precise, and domain agnostic interface for parsing
data encoded in the Crystallographic Information File (CIF) [@CIF] and Macromolecular
CIF (mmCIF) [@mmCIF] formats. Designed for programmatic analysis of crystalline systems,
`parsnip` offers a scriptable query API and a suite of convenient data retrieval methods
for automated lookups of structural information. Its minimal dependency set and
high-level interface allows for the rapid incorporation of `parsnip` into existing
libraries within the molecular simulation ecosystem [@Freud2020].

`parsnip`'s primary functionality lies in its ability to accurately reconstruct unit
cells from experimental data recorded with limited precision, often with only a few
decimal places. Through a combination of rational and floating-point arithmetic, we
achieve class-leading accuracy when reconstructing large and complex structures.
`parsnip`'s detailed processing of structural data yields better alignment with reported
space group and point group symmetries than existing tools [@ASE; @PyCIFRW; @GEMMI],
providing an ideal foundation for studies centered on material design.

`parsnip` supports a dictionary-like lookup format for both scalar and tabular data,
both of which can be expanded with Unix-style wildcards to simplify complex queries.
Convenient methods for parsing unit cell parameters, reconstructing particle positions,
and identifying site symmetry data are exposed to streamline common workflows in
materials data science. Parsnip's clear documentation of conventions and units
eliminates ambiguities common to interdisciplinary research. Using NumPy structured
arrays for data storage allows Python, C, and FORTRAN libraries to work with `parsnip`,
resulting in a stable, scalable dependency for scientific codebases in materials
research at the atomic, molecular, and colloidal scales [@Freud2020].

# Statement of Need

Materials scientists performing experimental and simulation research are fundamentally
exploring many of the same research questions. However, crystallographic software
designed for experimental data often does not scale well to automated workflows — a
particularly significant problem in interdisciplinary research where the building blocks
of crystal structures include atoms, macromolecules, and nanoparticles. While many
excellent libraries provide high-level interfaces and strict class hierarchies for
crystallographic data, the general nature of simulation science drives a need for
array-formatted storage that easily translates across simulation frameworks and system
length-scales. `parsnip` addresses this need by providing a simple, intuitive, and
well-documented software frontend that integrates tightly with existing standards for
molecular simulation and analysis. This marks a contrast in design between `parsnip` and
existing crystallography libraries like **ASE**, which provides wrapper types specific
to elemental systems and **PyCIFRW**, which lacks clear API documentation and lays out
data in a non-contiguous manner [@ASE; @PyCIFRW].

`parsnip` supports colloidal and mesoscale materials research in addition to the atomic
and protein datasets for which the CIF and mmCIF specifications were originally
designed. While the **gemmi** library also supports macromolecular data, `parsnip`
standardizes the data structures and API to ensure all inputs are handled in a
consistent, programmatic way [@GEMMI]. Rather than associating pure crystallographic
data with atomic symbols or valence states by default, `parsnip` provides only the
information required to reconstruct a particular structure unless otherwise queried.
This generality enables the application of decades of atomic and molecular research data
to novel studies of colloidal and soft matter crystallography. We also support
Unix-style wildcard queries, simplifying common lookup patterns like cell parameter
extraction and space group identification. Single-character wildcards enable
specification-compliant queries into heterogeneous databases of both CIF and mmCIF
files, further accelerating programmatic materials exploration. Although the Gemmi
library does support a similar style of wildcard through their `gemmi grep` command-line
tool, its use is limited to bash scripting and each wildcard query requires the file to
be re-parsed in its entirety [@GEMMI].

`parsnip` is designed and optimized for integration with larger materials science codes,
with a minimal dependency set and pure Python interface. This has facilitated its
incorporation into the **freud** analysis library, which uses `parsnip` to extract unit
cell data for use as reference structures in high-throughput simulation analysis
[@Freud2020]. Although **gemmi** also requires few dependencies, the compiled nature of
the package means that integration with other libraries may not be straightforward.
`parsnip` uses NumPy structured arrays to provide a stable data layout for
cross-language access without modifying the build system of downstream projects. As a
result, users get the benefits of copy-free data transfer to compiled languages without
any additional complexity.

# Acknowledgments

This research was supported by a Vannevar Bush Faculty Fellowship sponsored by the
Department of the Navy, Office of Naval Research under ONR award number
N00014-22-1-2821. Any opinions, findings, and conclusions or recommendations expressed
in this material are those of the authors and do not necessarily reflect the views of
the Office of Naval Research.

# References
