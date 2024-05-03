---
title: 'ASIMTools: A lightweight workflow manager for reproducible atomic simulations'
tags:
  - Python
  - atomic simulation
  - density functional theory
  - workflow
authors:
  - name: Mgcini Keith Phuthi
    orcid: 0000-0002-0982-8635
    equal-contrib: false
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Emil Annevelink
    orcid: 0000-0001-5035-7807
    equal-contrib: false
    affiliation: "2" # (Multiple affiliations must be quoted)
affiliations:
 - name: University of Michigan
   index: 1
 - name: Carnegie Mellon University
   index: 2
date: 20 October 2023
bibliography: paper.bib

# # Optional fields if submitting to a AAS journal too, see this blog post:
# # https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Atomic SIMulation Tools (`ASIMTools`) is a lightweight workflow and simulation manager for reproducible atomistic simulations that can be transferred across environments, DFT codes, interatomic potentials and structures implemented in Python for Unix based systems. By using in-built or user-defined asimmodules and utilities, users can run simulation recipes and automatically scale them on slurm based clusters or locally on their console. The core idea is to separate the dependence of the atomistic potential/calculator, the simulation environment and the simulation protocol thereby allowing the same simulation to be run with different calculators, structures or on different computers with just a change of one parameter in an input file after initial setup. This is extremely useful at a time when benchmarking Machine Learning Interactio Potentials has become a core part of computational materials science.Input and output yaml files follow a simple standard format, usually yaml, providing a simple interface that also acts as a record of the parameters used in a simulation without having to edit python scripts. The minimal set of requirements means any materials science codes can be incorporated into an ASIMTools workflow in a unified way.

# Statement of need

Atomic simulations are a key component of modern day materials science in both
academia and industry. However, simulation protocols and workflows used by
researchers are typically difficult to transfer to systems using different
inputs, codes and environments. It often involves rewriting entire scripts in
different languages to change from one type of atomistic potential or atomic
structure to another. This leads to poor reproducability and inefficient
transfer of code from one researcher to the next. In addition, there exists a
zoo of tools and packages for atomic simulation with more being developed every
day `[Walsh:2024]`. There is however no unifying framework that can encompass
all these tools without significant software development effort. Significant
effort should not be necessary because while the source of the fundamental
outputs of atomistic potentials such as energy, forces etc. may differ,
simulations built on these outputs should converge towards the most accurate
and computationally efficient. ASIMTools focuses on this last aspect by
introducing asimmodules which are simply Python functions that act as
simulation protocols which have no dependence on a specific atomistic potential
or computational environment or atomic structure. Through iteration and
community input, these simulation protocols will hopefully converge towards
best practice and ensure reproducibility of simulation results.

`ASIMTools` is for users interested in performing atomistic calculations on
UNIX-like operating systems and/or on slurm based High Performance Computing
clusters. By defining simulation protocols as functions in "asimmodules",
simulation protocols can be easily added to the library of provided asimmodules
and iterated on. This will allow the community to develop a robust set of
shareable simulation protocols. The flexibility of ASIMTools allows integration
of any kind of simulation tools such as pymatgen, LAMMPS etc. with examples
provided. With the asimmodules defined, users only need to provide a set of
inputs in the form of yaml files that define the parameters used for each
simulation and are therefore a record of these parameters. 

# State of the Field
There exist a number of popular workflow tools for atomistic simulations such
as Aiida `[@author:2001]`, Fireworks `[@author:2001]` and many more. These
tools provide frameworks for constructing complex workflows with different
underlying principles. Some managers enforce strict rules that ensure that data
obeys FAIR principles and emphasize data provenance and reproducibility. These
methods however tend to be fairly large packages with steep learning curves.
ASIMTools provides a simple interface as a starting point that can transform
any code into ASIMTools compatible code by simply wrapping it in a function
that returns a Python dictionary. Any such code can work in ASIMTools and with
a few extra steps, the protocol can be made to support an arbitrary calculator
and input structure.

In some workflow managers, such as Atomic Simulation Recipes `[@author:2001]`,
once workflows are built, it can often be difficult to quickly change and
iterate over key parameters such as the choice of atomistic calculator or
structure as they are intrinsically built into the code. This is particularly
challening in an age where machine learning models are becoming more popular.
Workflows involving machine learning interaction potentials tend to require the
ability to repeat the same calculations on different examples, using different
calculators on different hardware iteratively. This is where the value of
ASIMTools lies in contrast to more established workflows. ASIMTools is not
designed to replace the more powerful workflow managers but rather to
supplement them. This is achieved by providing unified inputs that can be
easily integrated into, for example, Aiida as Python functions/asimmodules
while also being a stand-alone lightweight workflow manager.

# Example
We present two examples of simulation protocols, many more can be found in the
ASIMTools documentation.

## Calculating the energy and forces of an atomic configuration
Most atomic simulations involve evaluations of energies, forces, dipoles etc.
of an atomic configuration. In Figure. \autoref{fig:singlepoint} we show how
the `singlepoint` asimmodule, provided in ASIMTools can be used and the input
files needed to run the asimmodule with arbitrary input structure, calculator
or environment.

![Schematic showing the connection between the modular input yaml files. The
sim_input.yaml is the main imput file which specifies the environment,
calculator (if used) and asimmodule to be
run.\label{fig:singlepoint}](figures/singlepoint.pdf){ width=100% }

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferredp
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Author Contribution Statement

Conceptualization by Keith Phuthi. Coding and development by Keith Phuthi and Emil Annevelink. Paper writing by Keith Phuthi and Emil Annevelink. Project management by all.
<!-- # Figures

Figures can be included like this:
![Schematic showing the connection between the modular input yaml files. The sim_input.yaml is the main imput file which specifies the environment, calculator (if used) and asimmodule to be run.\label{fig:singlepoint}](figures/singlepoint.pdf){ width=100% }
and referenced from text using \autoref{fig:example}. -->

# Acknowledgements

We acknowledge feedback from Kian Pu, Lance Kavalsky, Ziqi Wang and Hancheng Zhao.

# References
Example paper.bib file:

@article{Pearson:2017,
  	url = {http://adsabs.harvard.edu/abs/2017arXiv170304627P},
  	Archiveprefix = {arXiv},
  	Author = {{Pearson}, S. and {Price-Whelan}, A.~M. and {Johnston}, K.~V.},
  	Eprint = {1703.04627},
  	Journal = {ArXiv e-prints},
  	Keywords = {Astrophysics - Astrophysics of Galaxies},
  	Month = mar,
  	Title = {{Gaps in Globular Cluster Streams: Pal 5 and the Galactic Bar}},
  	Year = 2017
}

@book{Binney:2008,
  	url = {http://adsabs.harvard.edu/abs/2008gady.book.....B},
  	Author = {{Binney}, J. and {Tremaine}, S.},
  	Booktitle = {Galactic Dynamics: Second Edition, by James Binney and Scott Tremaine.~ISBN 978-0-691-13026-2 (HB).~Published by Princeton University Press, Princeton, NJ USA, 2008.},
  	Publisher = {Princeton University Press},
  	Title = {{Galactic Dynamics: Second Edition}},
  	Year = 2008
}

@article{Walsh:2024,
	Title = {Open computational materials science},
	Volume = {23},
	copyright = {2024 Springer Nature Limited},
	issn = {1476-4660},
	url = {https://www.nature.com/articles/s41563-023-01699-7},
	doi = {10.1038/s41563-023-01699-7},
	Abstract = {The materials modelling community is emerging as a champion for reproducible and reusable science. Aron Walsh discusses how FAIR databases, collaborative codes and transparent workflows are advancing this movement.},
	Language = {en},
	Number = {1},
	urldate = {2024-01-11},
	Journal = {Nature Materials},
	Author = {Walsh, Aron},
	Month = jan,
	Year = {2024},
	Note = {Number: 1
  Publisher: Nature Publishing Group},
	Keywords = {Research data, Theory and computation},
	Pages = {16--17},
}

