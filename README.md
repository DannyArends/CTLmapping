R package CTL
================
CTL provides an implementation in R for the novel Correlated Trait Locus (CTL) mapping
mapping methodology. CTL mapping is a novel approach to detect genetic regulation of 
phenotypes in recombinant inbred line populations (RIL). It is a method which complements
classical QTL analysis, providing additional insights overlooked by the classical QTL 
approach. 

Differences in correlations between traits within an inbred population are determined 
at each genetic marker. Phenotypes are assigned to genotype groups and a single phenotype 
is used to scan all other phenotypes for a loss in or a gain of correlation. The likelihood 
profiles (similar to QTL profiles) of this 'loss of correlation' measurement show a very 
high degree of overlap with classical QTL profiles, BUT additional information is easily 
extracted from phenotype x phenotype interaction plots. With the right dataset (ideally a 
combination of: classical phenotypes, protein abundance and gene expression) CTL shows 
the genetic wiring of the classical phenotypes and identify key players in the genetic / 
protein network underlying QTL and CTL.

Installation
------------
Prepare your environment by following these steps:

- Download and Install the R environment from [www.r-project.org](http://www.r-project.org/ "www.r-project.org")

Then install into R by using (from a terminal / command line):

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ R CMD INSTALL CTLmapping                               # Install the package

Optionally you can install the pre-build packages by downloading the appropriate 
package for your operating system. 

Starting
--------
Load the library in the R interface by the following command (in R):

```R
    library(ctl)                            # Load the library
    ?ctl                                    # Show the help
```

Examples
========
Scan your data

```R
    library(ctl)
    data(multitrait)
    multitrait = fill.geno(multitrait)      # Fill missing genotype values
    ?CTLscan                                # Show the CTLscan function help
    ctl_result <- CTLscan.cross(multitrait)
```

Plot a single phenotype, the profile is comparable to the QTL profile. However using 
CTL mapping we know which phenotypes are differentially correlated underneath the peak.
This additional information adds to the already known QTL information.

```R
    plot(ctl_result, pheno.col=12)
```

Create an image of the phenotypes to marker relation strength, this matrix is 'comparable' 
to a heat map of QTL scans on many phenotypes, the underlying model assumptions are different 
from QTL mapping but comparable, thus the output is not shockingly different from QTL mapping.

```R
    r1 <- image(ctl_result,against="markers")
```

Create an image of the phenotypes to phenotypes relation strength, this is the additional 
information matrix, which is not available in classical QTL mapping.

```R
    r2 <- image(ctl_result,against="phenotypes")
```

Reconstruct the network and write two sif files. One sif file contains the full network, the other 
holds the edge summary network.

```R
    $ > CTLnetwork(ctl_result)
```

We can use Cytoscape to visualize the created network (available from [www.cytoscape.org](http://www.cytoscape.org// "www.cytoscape.org") )

Contributing and TODO
---------------------
Want to contribute? Great!

See inst/TODO.txt

Just post comments on code / commits.
Or be a maintainer, and adopt a function

Danny Arends

Disclaimer
----------
Copyright (c) 2010 Danny Arends
