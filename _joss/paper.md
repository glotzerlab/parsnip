---
title: "parsnip: Streamlined Crystallographic Data Parsing for Simulation Science"
tags:
  - Python
  - crystallography
  - materials science
authors:
  - name: Jenna Bradley
    orcid: 0009-0007-2443-2982
    affiliation: 1
  - name: Sharon C. Glotzer
    orcid: 0000-0002-7197-0085
    affiliation: "[1, 2]"
affiliations:
 - name: Materials Science and Engineering, University of Michigan, United States
   index: 1
   ror: 00jmfr291
 - name: Department of Chemical Engineering, University of Michigan
   index: 2
date: 2 January 2026
bibliography: paper.bib
---

# Summary

`parsnip` provides a lightweight, precise, and domain agnostic interface for parsing
data encoded in the Crystallographic Information File (CIF) [@CIF] and Macromolecular
CIF (mmCIF) [@mmCIF] formats. Designed for programmatic analysis of crystalline systems,
`parsnip` offers a scriptable API and a suite of convenient data retrieval methods for
automated lookups of structural information. Its minimal dependency set and high-level
interface allows for the rapid incorporation of `parsnip` into existing libraries within
the molecular simulation ecosystem [@Freud2020].

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
materials data science. `parsnip`'s clear documentation of conventions and units
eliminates ambiguities common to interdisciplinary research. Our use of NumPy structured
arrays for data storage allows Python, C, and FORTRAN libraries to work with `parsnip`,
resulting in a stable, scalable dependency for scientific codebases in materials
research at the atomic, molecular, and colloidal scales.

# Statement of Need

<!-- A section that clearly illustrates the research purpose of the software and places it in the context of related work. This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work. -->

More than thirty years of material data of is encoded in the CIF file format, with
terabytes of elemental and protein structures freely available to researchers
[@MATERIALSPROJECT; @PDB]. While early CIF parsers were predominantly written in C and
Fortran, the growth of Python opened new opportunities for simple, scriptable access to
crystallographic data. **PyCifRW** is one of the earliest such tools, offering a
complete and specification complicant parser for the CIF format [@PyCIFRW]. Soon
thereafter, **ASE**, a tool designed to initialize atomistic simulations, added support
for CIF files as an alternative to NetCDF inputs for DFT simulations [@ASE_ORIG; @ASE].
This marked a transition from pure IO libraries to mixed IO and analysis tools, with
both **ASE** and later projects like **pymatgen** including code to characterize and
manipulate structures [@pymatgen]. One of the most recent CIF libraries, **gemmi**, also
follows this trend, although significant amounts of code are committed to quickly and
accurately parsing the CIF grammar.

While all of these tools provide excellent interfaces for researchers working with
atomic materials, the structure and typing of alternative libraries is ill-suited to
interdisciplinary research, where the building blocks of crystal structures include
atoms, macromolecules, and nanoparticles. The explosion of simulation research in
superatomic systems has driven a need for scalable, array-formatted crystallographic
data that easily translates across simulation frameworks and system length-scales.
`parsnip` addresses this need by providing a simple, intuitive, and well-documented
software frontend that integrates tightly with existing standards for molecular
simulation and analysis [@HOOMD; @LAMMPS].

# State of the Field

<!-- (~200 words): A description of how this software compares to other commonly-used packages in the research area. If related tools exist, provide a clear "build vs. contribute" justification explaining your unique scholarly contribution and why existing alternatives are insufficient. -->

`parsnip` supports colloidal and mesoscale materials research in addition to the atomic
and protein datasets for which the CIF and mmCIF specifications were originally
designed. This is a fundamental change to the scope of the library, and necessitates
novel tooling to support the greater diversity of systems. For example, while both
`parsnip` and the **gemmi** [@GEMMI] library support macromolecular data encoded in the
mmCIF format, `parsnip` standardizes its API to ensure all inputs are handled in a
consistent, programmatic way. In contrast, **gemmi** separates the parsing API and data
structures between the two, improving performance at the expense of generality.

Existing crystallography libraries like **ASE** [@ASE] and **pymatgen** [@pymatgen]
encode atomic information into the data structures and types of parsed information,
requiring postprocessing for studies where that information is unnecessary or incorrect.
Rather than associating pure crystallographic data with atomic symbols or valence states
by default, `parsnip` provides only the information required to reconstruct a particular
structure unless otherwise queried. This approach offers benefits for users who do
require atomic information, as we are able to provide access to arbirary data associated
with the basis positions, rather than fixed keys as required by other tools.

# Software Design

<!-- Explain the trade-offs you weighed, the design/architecture you settled on, and why it matters. If other packages exist in this space, you must include an explicit "build vs. contribute" justification explaining why you created new software rather than contributing to existing projects. -->

`parsnip`'s design enables studies of crystallinity in colloidal matter, a feature set
that is not met by any available tools in the field. By separating units and particle
data from pure structural information, we are able to provide a general interface for
the study of ordered matter across scales. We also abstract away some portions of the
CIF syntax to simplify the API -- most notable, we aggregate across `data_` blocks to
ensure all relevant information is accessible through a uniform grammar of queries. This
choice enables a user interface that more closely resembles other structural data
formats like *XYZ*, *MOL*, and *VTP* [@molIUPAC; @vtkBook].

Our unique approach to the CIF specification extends to the design of our parser as
well. While most existing tools in the space use parser generators based on the IUCR's
formal grammar, we identified a non-neglible subset of CIF files that break the formal
specification but nevertheless contain useful data. To overcome this, `parsnip` does not
validate the entire syntax tree of the CIF grammar: rather, we eagerly consume nodes
near the leaves of the tree that appear to contain data. This is a departure from the
standard "validating" parser strategy, but it enables fast and robust data extraction
without a significant increase in code complexity.

While this parsing technique alone provides significant accuracy benefits, there are
still many files that cannot be accurately reconstructed by other tools. Standard,
"symmetrized" CIF data requires the application of symmetry operations to reconstruct a
lattice. These operations have the form of the form $(x \pm \mathbb{Q}) \mod 1$, where
$x$ is experimentally determined. While Wyckoff positions can have arbitrary real
values, the data stored in CIF files is necessarily finite. For this reason, the actual
set of valid, parseable positions is the group $\mathbb{Q} \mod 1$. Rather than
evaluating expressions in floating point arithmetic like other CIF libraries, `parsnip`
actually evaluates unit cell positions in the correct rational form. We then convert to
back to floating point values for a tolerance-based deduplication check, which catches
edge cases in recorded data where values are not rounded consistently (e.g.
$\left(1/3, 2/3\right)=\left(0.3333, 0.6666\right)$).

# Research Impact Statement

<!-- Provide concise evidence of either realized impact (e.g., external use, integrations, enabled results) or credible near-term significance (novel capability with benchmarks, reproducible materials, and community-readiness signals such as documentation, tests, license, releases, and contribution process). -->

`parsnip` is designed and optimized for integration with larger materials science codes,
with a minimal dependency set and pure Python interface. This has facilitated its
incorporation into the **freud** analysis library, which uses `parsnip` to build
reference structures for high-throughput simulation analysis [@Freud2020]. When studying
complex crystals, standard characterization techniques often fail to uniquely resolve a
structure of interest. Simplifying access to a massive array of crystal structures
stored in databases like the Crystallography Open Database (COD) [@COD] enables the
construction of reference datasets for both classical and machine-learned
characterization techniques.

`parsnip` also supports Unix-style wildcard queries, simplifying common lookup patterns
like cell parameter extraction and space group identification. Single-character
wildcards enable specification-compliant queries into heterogeneous databases of both
CIF and mmCIF files, further accelerating programmatic materials exploration. Although
the **gemmi** library does support a similar style of wildcard through their
`gemmi grep` command-line tool, its use is limited to bash scripting and each wildcard
query requires the file to be re-parsed in its entirety [@GEMMI].

Tests against 10,099 CIF files from the Crystallography Open Database (COD) shows we are
able to correctly extract 95.4% of structures, more than any other library we could
find. Table \ref{accuracyCOD} shows `parsnip`'s excellent performance compared to its
contemporaries: **parnsnip**'s symbolic parsing mode is the most accurate of all tested
CIF libraries, and is able to correctly reconstruct more files than the next best
alternative, **ASE**. We note that our results use a single, fixed parsing precision for
all 10,099 files: however, as discussed `in parsnip`'s documentation, tailoring the
parse precision to match the precision of the data in the file yields even better
results.

: Comparison of unit-cell reconstruction accuracy for 10099 CIF files from the
COD.\label{accuracyCOD}

| Library                   | Incorrect Structure | Failed Parses | Total Errors |
| ------------------------- | :-----------------: | :-----------: | :----------: |
| `parsnip (sympy, prec=3)` |       **84**        |      385      |   **469**    |
| `parsnip (prec=3)`        |         203         |      385      |     588      |
| ASE                       |         125         |      810      |     935      |
| gemmi                     |        2343         |     **0**     |     2343     |

<!-- `parsnip` is also extremely well documented, with tutorials building from a simple API -->
<!-- introduction all the way to a detailed guide on using the code with the `HOOMD-Blue` and -->
<!-- `LAMMPS` simulation engines. We also illustrate the use of `parsnip`'s symbolic parsing -->
<!-- functionality, with examples on how to ensure the correctness of structures built from -->
<!-- low-quality data. This is a notable improvement from older libraries like `pycodcif` -->
<!-- [@CITE] and `PyCIFRW` [@PyCIFRW], which have limited Python API documentation and few -->
<!-- examples. -->

# Acknowledgments

This research was supported by a Vannevar Bush Faculty Fellowship sponsored by the
Department of the Navy, Office of Naval Research under ONR award number
N00014-22-1-2821.

# AI Usage Disclosure

GitHub Copilot’s automated PR review was evaluated for use with this project, but was
disabled due to the unhelpful nature of the output. No generative AI tools were used in
the development of this software, the writing of this manuscript, or the preparation of
supporting materials.

# References
