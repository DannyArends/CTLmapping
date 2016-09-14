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

CTL mapping (CTL) is a novel method to associate correlation differences observed between phenotypes to genomic locations. 
CTL mapping can be used in the analysis of experimental and outbred crosses, and can be used to study human populations.

CTL allows geneticists, biologists and animal breeders to analyze correlation difference on phenotypes covering the whole 
biomolecular spectrum. From 'classic' phenotypes such as yield and disease suscepibility to high-throughput experimental 
data such as microarrays, RNA-seq and protein abundance measurements.

The CTL software is provided as a free and open source (FOSS) package for the R Project for Statistical Computing [@R:2005]. The 
core algorithm is written in C allowing it to be deployed anywhere, but it also provides easy integration of the algorithm 
into other languages. As a proof of concept we provide bindings for the D 2.0 language.

Data structures of the CTL mapping package have been harmonized with the popular R/qtl package [@Arends:2010], allowing 
users to quickly and effiently re-analyse previous (R/)QTL experiments. Additional advantages of close integration with 
R/qtl are that CTL users can use any input formats supported by R/qtl, and additionally have access to all plot and 
helper functions provided by R/qtl.

Results from CTL mapping can be visualized by several plotting routines provided by the package, however results can 
also be exported to external tools (such as cytoscape) for visualization.

# Example datasets

CTL mapping comes with several example datasets (in Rdata format) for the user to explore:

- 301 gene expression traits measured on 109 Saccharomyces cerevisia [@Brem:2002]
- 9 Metabolite expression traits measured on 403 Arabidopsis Thaliana [@Churchill:2012]
- 24 Metabolite expression traits measured on 162 Arabidopsis Thaliana [@Keurentjes:2002]

-![CTL X3.Methylthiopropyl](CTL.png)

# Future work

Future work includes research into better handling of RAM by the CTL algorithm. 
The authors also aim to provide additional interactive visualization (such as D3 interactive graphics). 

# References