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

# Summary

Parsnip provides a lightweight and well-documented interface for reading
Crystallographic Information Files (CIFs).

# Statement of need

Materials scientists performing experimental or simulation research are fundementally
exploring many of the same research questions. However, the two parties differ in the
preferred structure and analysis techniques, requi <!--TODO -->. While many excellent
libraries provide high-level interfaces and strongly typed class hierarchies for
crystallographic data, the general nature of simulation science drives a need for
array-formatted data with intuitive memory layouts that easily translate between
simulation frameworks.

`parsnip` is targeted at colloidal and mesoscale materials research, rather than the
atomic and protein datasets that the CIF and mmCIF specifications were originally
targeted at. Rather than associating pure crystallographic data with atomic symbols or
valence states by default, `parsnip` provides only the information required to
reconstruct a particular structure unless otherwise queried. This generality allows
users unfamiliar with the CIF specification to process crystallographic information
files, as dedicated lookup methods automatically locate the correct queries that return
unit cell and basis site information. We also support Unix-style wildcard queries, a
feature that is not included in the Python APIs of any similar package. This feature
allows for rapid access to "slices" of the CIF keyword specification, simplifying common
complex lookups to much simpler queries.

`parsnip` is designed and optimized for use as a dependency in larger materials science
codes, requiring only Numpy as a dependency. For this reason, it has already been
incorporated into the **Freud** analysis library, which uses `parsnip` to extract unit
cell data for use as reference structures in high-throughput simulation analysis.
Although other CIF libraries like **Gemmi** have minimal dependency sets, the compiled
nature of the underlying library means that integration with other libraries may be
nontrivial. `parsnip` uses Numpy structured arrays to provide a stable, in memory data
layout for cross-language access without modifying the build system of downstream
projects. As a result, users get the benefits of copyless data transfer to compiled
languages without any additional complexity.

<!-- TODO: symbolic parsing -->

# Examples

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations \begin{equation}\label{eq:fourier} \hat
f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx \end{equation} and refer to
\autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a
preferred citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:

- `@author:2001` -> "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png) and referenced from text
using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

# References
