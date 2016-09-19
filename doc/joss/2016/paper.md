---
title: 'CTL: Discovering genetic loci associated with correlation differences'
tags:
  - bioinformatics
  - genetics
  - genomics
authors:
  - name: Danny Arends
    orcid: 0000-0001-8738-0162
    affiliation: Humboldt University, Berlin, Germany, University of Groningen, The Netherlands
  - name: Yang Li
    orcid:
    affiliation: University of Groningen, The Netherlands, University Medical Center Groningen, The Netherlands
  - name: Gudrun A. Brockmann
    orcid:
    affiliation: Humboldt-Universität zu Berlin, Germany
  - name: Ritsert C. Jansen
    orcid:
    affiliation: University of Groningen, The Netherlands
  - name: Robert W. Williams
    orcid: 0000-0001-8924-4447
    affiliation: University of Tennessee Health Science Center, USA
  - name: Pjotr Prins
    orcid: 0000-0002-8021-9162
    affiliation: University Medical Center Utrecht, The Netherlands, University of Tennessee Health Science Center, USA
date: 29 May 2016
bibliography: paper.bib
---

# Summary

CTLmapping is an implementation of the Correlation Trait Loci (CTL)
algorithm first presented in [@Arends:thesis_chapter]. CTLmapping
allows geneticists to analyze correlation
difference between phenotypes.  CTLmapping is complementary to the
proven quantitative trait locus (QTL) mapping method which correlates
observed phenotype against genotype. CTL mapping
associates correlation differences observed *between* phenotypes,
subject to the genotype. In other words, QTL mapping treats phenotypes
independently while CTL mapping connects phenotypes. CTL generally
show very similar profiles to QTL, but get interesting when they
differ (see figure 1).

Because CTL connect phenotypes CTLmapping provides a mechanism for
inference and discovering causality [@Arends:thesis_chapter].  This is
particularly of interest when phenotype correlations change with
conditions, for example in pathways with highly correlated gene
expression patterns (see figure 1).  CTLmapping differs from set test
methods, e.g., [@Wang:2010] in that CTLmapping does not require prior information on sets
(e.g., pathways).

-![Figure 1](Fig1.png)
Figure 1: CTL profile and QTL profiles. (a) shows a very similar
profile, (b) shows a CTL without a QTL and (m) shows a QTL without a
CTL.

CTL analysis can be performed on phenotypes obtained from the whole
biomolecular spectrum. From 'classic' phenotypes, such as yield and
disease suscepibility, to high-throughput experimental data, such as
microarrays, RNA-seq and/or protein abundance measurements.  This is
especially useful in combined datasets, e.g. a combination of:
classical phenotypes, protein abundance and gene expression (see
figure 2).

-![Figure 2](Fig2.png)
Figure 2: CTL show the genetic wiring of classical phenotypes and
identify key players in the genetic / protein network underlying
classical phenotypes using QTL and CTL information.

CTLmapping can be applied in model organism experimental and outbred
crosses, such as mouse and the plant *Arabidopsis thaliana* (see
example datasets below), as well as in natural populations, such as
human. For statistical power the more individuals the better, but as a
rule of thumb it is about the same as for QTL, i.e., about 100
individuals for a recombinant inbred line CTL, and 1,000 individuals
for GWAS-style CTL [@Arends:thesis_chapter].

The CTLmapping software is provided as a free and open source (FOSS)
package for the R Project for Statistical Computing [@R:2005].
Data structures of the CTL mapping R package have been harmonized with
the popular R/qtl package [@Arends:2010], allowing users to quickly
and efficiently re-analyse previous (R/)QTL experiments. Additional
advantages of close integration with R/qtl are the many input formats
supported by R/qtl, and access to all plot and helper functions
provided by R/qtl.

The core CTLmapping algorithm is written in standalone C making it
easy to integrate the CTL mapping algorithm into other languages that
allow calling C functions. As a proof of concept the CTL repository
provides bindings for the D language.

CTL has been integrated into GeneNetwork (GN), a FOSS framework for
web-based genetics that can be deployed anywhere [@Sloan:2016]. This
allows results from CTL mapping to be interactively explored using the
GeneNetwork web interface. Additionally results from CTL mapping can
be visualized by plotting routines provided by the R package and
results can be exported to external tools (such as Cytoscape
[@Cytoscape:2003]) for visualization and interactive
exploration.

# Example datasets

CTL mapping comes with several example datasets (in Rdata format) for
the user to explore:

- 301 gene expression traits measured on 109 *Saccharomyces cerevisia* [@Brem:2002]
- 9 Metabolite expression traits measured on 403 *Arabidopsis Thaliana* [@Churchill:2012]
- 24 Metabolite expression traits measured on 162 *Arabidopsis Thaliana* [@Keurentjes:2002]

(instructions can be found in the README).

# Future work

CTL is computationally very intensive, phenotypes O(n^2), both in
terms of RAM use and CPU.  Future work includes research into
improving the CTL algorithm for large scale correlations and inference
using GPU/supercomputing.  We are also working on adding explorative
interactive visualization (such as Cytoscape and D3 interactive
graphics).

# References
