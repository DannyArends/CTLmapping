## Correlated Trait Locus (CTL) mapping

In this repository implementations of the novel Correlated Trait Locus (CTL) 
mapping algorithm. Provided as package for the statistical language R and in 
the D 2.0 programming language.

CTL mapping is a novel approach to detect genetic regulation of phenotypes in 
natural and experimental populations. It is a method which complements classical 
QTL analysis, providing additional insights overlooked by the classical QTL 
approach.

### Algorithm

Differences in correlations between traits within an inbred population are 
determined at each genetic marker. Phenotypes are assigned to genotype groups 
and a single phenotype is used to scan all other phenotypes for a loss or gain 
of correlation. The likelihood profiles (~ QTL profiles) of this 'loss of 
correlation' measurement shows a very high degree of overlap with classical 
QTL profiles. However additional information is available from phenotype x 
phenotype interactions. With the right dataset (ideally a combination of: 
classical phenotypes, protein abundance and gene expression) CTL shows the 
genetic wiring of the classical phenotypes and identify key players in the 
genetic / protein network underlying QTL and CTL.

### Download the package
The first and quickest way to start mapping CTLs, is to install the package 
directly into R from Github using the devtools package:

```R
# install.packages("devtools")
library(devtools)
install_github("CTLmapping", "DannyArends", subdir="Rctl")
```

The second option to is to clone the package from Github. First prepare your 
environment by download and 'moving' to the folder:

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping                                          # Goto the folder

### Use the R library

Prepare your environment by download and installing the R environment from 
[www.r-project.org](http://www.r-project.org/ "www.r-project.org"). Then 
download CTLmapping and install into R by using (from a terminal / command 
line):

    $ R CMD INSTALL Rctl                                     # Install the package

or use the 'installR' makefile target:

    $ make installR                                          # Install into R

Infuture packages will be on CRAN, and optionally installing a pre-build packages 
by downloading the appropriate package for your R version and operating system. 

### Compile the standalone executable

#### (C version)

Just run 'make' from a terminal / command line:

    $ make versionC                                          # Compile the executable
    $ make static                                            # Compile the static library
    $ make shared                                            # Compile the shared library

#### (D 2.0 version)

Prepare your environment by download and installing the DMD 2.0 compiler from 
[www.d-programming-language.org](http://www.d-programming-language.org 
"www.d-programming-language.org"). Run 'make' from a terminal / command line:

    $ make versionD                                          # Compile the executable

Optionally you can use a provided binary by downloading the approriate one for your 
operating system.

### Starting in R

Load the library in the R interface by the following command (in R):

```R
library(ctl)                            # Load the library
?ctl                                    # Show the help
```

### Examples

Scan your data

```R
library(ctl)
data(multitrait)
ctlres = CTLscan.cross(multitrait)
```

Plot a single phenotype, the profile is comparable to the QTL profile. However using 
CTL mapping we know which phenotypes are differentially correlated underneath the peak.
This additional information adds to the already known QTL information.

```R
plot(ctlres, pheno.col=12)
```

Create an image of the phenotypes to marker relation strength, this matrix is 'comparable' 
to a heat map of QTL scans on many phenotypes, the underlying model assumptions are different 
from QTL mapping but comparable, thus the output is not shockingly different from QTL mapping.

```R
r1 = image(ctlres,against="markers")
```

Create an image of the phenotypes to phenotypes relation strength, this is the additional 
information matrix, which is not available in classical QTL mapping.

```R
r2 = image(ctlres,against="phenotypes")
```

Reconstruct the network and write two sif files. One sif file contains the full network, the other 
holds the edge summary network.

```R
CTLnetwork(ctlres)
```

We can use Cytoscape to visualize the created network (available from [www.cytoscape.org](http://www.cytoscape.org// "www.cytoscape.org") )

### Contributing and TODO

Want to contribute? Great!

See Rctl/inst/TODO.txt or submit a Github issue or pull request

Its also possible to just post comments on code / commits.
Or be a maintainer, and adopt a function

Danny Arends

### Disclaimer

Copyright (c) 2010-2013 GBIC: Danny Arends, Ritsert C Jansen, Pjotr Prins

