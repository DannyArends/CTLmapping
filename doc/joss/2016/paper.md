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
    affiliation: Humboldt-Universit�t zu Berlin, Germany
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

CTL mapping is a novel method that associates correlation differences observed between phenotypes to genomic loci, 
These loci, called Correlated Trait Loci (CTL) are regions in the genome at which different genotypes show a significant 
difference in correlation. CTL mapping can be used in the analysis of experimental and outbred crosses, as well as 
human populations. CTL allows geneticists, biologists and animal breeders to analyze correlation difference between 
phenotypes, and exploit these difference in correlation in selective breeding to break adverse phenotype to phenotype 
correlations. CTL analysis can be performed on phenotypes obtained from the whole biomolecular spectrum. From 'classic' 
phenotypes such as yield and disease suscepibility to high-throughput experimental data such as microarrays, RNA-seq 
and/or protein abundance measurements.

The CTL software is provided as a free and open source (FOSS) package for the R Project for Statistical Computing 
[@R:2005]. The core algorithm is written in C allowing it to be deployed anywhere, but it also provides easy 
integration of the CTL mapping algorithm into other languages that allow calling C functions. As a proof of concept 
the CTL repository provide bindings for the D 2.0 language. 

Data structures of the CTL mapping R package have been harmonized with the popular R/qtl package [@Arends:2010], allowing 
users to quickly and effiently re-analyse previous (R/)QTL experiments. Additional advantages of close integration with 
R/qtl are the many input formats supported by R/qtl, and access to all plot and helper functions provided by R/qtl.

To increase the reach of our method, we have integrated CTL into GeneNetwork (GN) a FOSS framework for web-based 
genetics that can be deployed anywhere [@Sloan:2016]. This allows results from CTL mapping to be interactively 
explored using the familiar GeneNetwork web interface. Additionally results from CTL mapping can be visualized 
by plotting routines provided by the R package, however results can also be exported to external tools (such as 
Cytoscape [@Cytoscape:2003]) for visualization and interactive exploration.

# Example datasets

CTL mapping comes with several example datasets (in Rdata format) for the user to explore:

- 301 gene expression traits measured on 109 Saccharomyces cerevisia [@Brem:2002]
- 9 Metabolite expression traits measured on 403 Arabidopsis Thaliana [@Churchill:2012]
- 24 Metabolite expression traits measured on 162 Arabidopsis Thaliana [@Keurentjes:2002]

-![QTL profile (Positive Y-axis) and CTL associations (Negative Y-axis) for X3.Methylthiopropyl](CTL.png)

# Future work

Future work includes research into better handling of RAM by the CTL algorithm. 
The authors also aim to provide additional interactive visualization (such as D3 interactive graphics). 

# References