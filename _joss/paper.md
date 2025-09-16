---
title: "Parsnip: Streamlined Crystallographic Data Parsing for Simulation Science"
tags:
  - Python
  - crystallography
  - materials science
authors:
  - name: Jenna Bradley
    orcid: 0009-0007-2443-2982
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Materials Science and Engineering, University of Michigan, United States
   index: 1
   ror: 00jmfr291
date: XX November 2025
bibliography: paper.bib
---

<!-- Stronger than lightweight/well documented ->: accurate! (rational arithm.) (leightweight good. meat should be fn), correct? -->
<!-- empasize: still works for atomistic systems, but works well for non-atomistic systems .-->

<!-- easily queryable? intuitive?-->

# Summary

`parsnip` provides a lightweight, highly precise, and domain agnostic interface for
parsing material data encoded in the Crystallographic Information File (CIF) and
Macromolecular CIF (mmCIF) formats. Designed for programmatic analyses of crystalline
systems, `parsnip` offers a scriptable query API and a suite of convenient data
retrieval methods for automated lookups of structural information. Its minimal
dependency set and high-level interface has allowed for the rapid incorporation of
`parsnip` into several existing libraries in the molecular simulation ecosystem.

`parsnip`'s primary functionality lies in its ability to accurately reconstruct unit
cells from only a few decimal places of recorded experimental data. Through a
combination of decimal and floating-point arithmetic, we achieve class-leading accuracy
in reconstructing large and complex structures. The detailed processing of structural
data yields better alignment with reported space group and point group symmetries than
existing tools, providing an ideal foundation for studies centered on material design.

<!-- TODO: explain list of features, with code block example. one more paragraph-->

# Statement of Need

Materials scientists performing experimental and simulation research are fundamentally
exploring many of the same research questions. However, the two parties benefit from
specialized software tailored to the needs of research techniques. While many excellent
libraries provide high-level interfaces and strict class hierarchies for
crystallographic data, the general nature of simulation science drives a need for
array-formatted data with intuitive memory layouts that easily translate between
simulation frameworks. This shift in design focus provides a simple, intuitive software
frontend that integrates tightly with existing standards for molecular simulation and
analysis. This marks a contrast in design between `parsnip` and existing crystallography
libraries like **ASE**, which provides wrapper types specific to elemental systems and
**PyCIFRW**, which lays out data in a noncontiguous manner [@Larsen2017; @Hester2006].

`parsnip` targets colloidal and mesoscale materials research in addition to the atomic
and protein datasets that the CIF and mmCIF specifications were originally designed for.
Rather than associating pure crystallographic data with atomic symbols or valence states
by default, `parsnip` provides only the information required to reconstruct a particular
structure unless otherwise queried. This generality enables the application of decades
of materials research data to novel explorations of colloidal and soft matter
crystallography. We also support Unix-style wildcard queries, simplifying common lookup
patterns like cell parameter extraction (6 keys) and space group identification
(variable 1-5 keys) to general, reusable expressions. Single-character wildcards enable
specification-compliant queries into heterogeneous databases of both CIF and mmCIF
files, further accelerating programmatic materials exploration. Although the Gemmi
library does support a similar style of wildcard in their `gemmi grep` command-line
tool, it's use is limited to bash scripting and each wildcard query requires the file to
be re-parsed in its entirety [@Wojdyr:2022].

`parsnip` is designed and optimized for use as a dependency in larger materials science
codes, requiring only Numpy as a dependency. For this reason, it has already been
incorporated into the **Freud** analysis library, which uses `parsnip` to extract unit
cell data for use as reference structures in high-throughput simulation analysis.
Although other CIF libraries like **Gemmi** have minimal dependency sets, the compiled
nature of the underlying library means that integration with other libraries may be
nontrivial. `parsnip` uses Numpy structured arrays to provide a stable data layout for
cross-language access without modifying the build system of downstream projects. As a
result, users get the benefits of copyless data transfer to compiled languages without
any additional complexity.

<!-- TODO: show example code block as figure? -->
<!-- TODO: mention symbolic parsing : 'a primary challenge...'-->

<!-- # Citations -->

<!-- Citations to entries in paper.bib should be in -->
<!-- [rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html) -->
<!-- format. -->

<!-- If you want to cite a software repository URL (e.g. something on GitHub without a -->
<!-- preferred citation) then you can do it with the example BibTeX entry below for @fidgit. -->

<!-- For a quick reference, the following citation commands can be used: -->

<!-- - `@author:2001` -> "Author et al. (2001)" -->
<!-- - `[@author:2001]` -> "(Author et al., 2001)" -->
<!-- - `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)" -->

<!-- # Figures -->

<!-- Figures can be included like this: -->
<!-- ![Caption for example figure.\label{fig:example}](figure.png) and referenced from text -->
<!-- using \autoref{fig:example}. -->

<!-- Figure sizes can be customized by adding an optional second parameter: -->
<!-- ![Caption for example figure.](figure.png){ width=20% } -->

<!-- # Acknowledgements -->

# References
