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
QTL profiles. However additional information is available from the phenotype x 
phenotype interactions. With the right dataset (ideally a combination of: 
classical phenotypes, protein abundance and gene expression) CTL shows the 
genetic wiring of the classical phenotypes and identify key players in the 
genetic / protein network underlying classical phenotypes using QTL and CTL 
information.

### Installing the R package
The quickest way to start mapping CTLs, is to install the package directly into 
R from Github using the devtools package:

```
# install.packages("devtools")
library(devtools)
install_github("CTLmapping", "DannyArends", subdir="Rctl")
```

After installing [learn the R commands](#starting-in-r "Starting in R")

### Download the software (R, C and D)

The second option to is to clone the package from Github, this will give you 
access to the R package but also provide the command line tools. First prepare 
your environment by download and 'moving' to the folder:

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

Plans are to put the package on CRAN, but this has not happend yet. 

### Compile the standalone executable

Optimized versions of the software are also available for more high throughput data.
Here we explain how to build a standalone executable version using C

#### (C version)

Just run 'make' from a terminal / command line:

    $ make versionC                                          # Compile the executable
    $ make static                                            # Compile the static library
    $ make shared                                            # Compile the shared library

C code can also be compiled also into static or dynamic link libraries, when using the 
appropriate makefile commands. 

#### (D 2.0 version)

Prepare your environment by download and installing the DMD 2.0 compiler from 
[www.d-programming-language.org](http://www.d-programming-language.org 
"www.d-programming-language.org"). Run 'make' from a terminal / command line:

    $ make versionD                                          # Compile the executable

Optionally you can use a provided binary by downloading the approriate one for your 
operating system.

### Starting in R

Load the library in the R interface by the following command (in R):

```
library(ctl)                            # Load the library
?ctl                                    # Show the help
```

### Examples

If you have data prepared for (R/qtl)[http://www.rqtl.org/ "R/qtl website"], 
you can use the CTLscan.cross() function to directly scan your R/qtl cross 
object.

Scan some example data (in R/qtl format) using the CTLscan.cross function:

```
library(ctl)
data(multitrait)
ctlres = CTLscan.cross(multitrait)
```

There is also a CTLscan() function, this function takes plain old genotype and 
phenotype matrices as input, and is can be called in the same way.

```
library(ctl)
data(ath.metabolites)
ctlres = CTLscan(ath.metab$genotypes, ath.metab$phenotypes)
```

Plot a single phenotype, the profile is comparable to the QTL profile. However using 
CTL mapping we know which phenotypes are differentially correlated underneath the peak.
This additional information adds to the already known QTL information.

```
plot(ctlres, pheno.col=12)
```

Create an image of the phenotypes to marker relation strength, this matrix is 'comparable' 
to a heat map of QTL scans on many phenotypes, the underlying model assumptions are different 
from QTL mapping but comparable, thus the output is not shockingly different from QTL mapping.

```
r1 = image(ctlres, against="markers")
```

Create an image of the phenotypes to phenotypes relation strength, this is the additional 
information matrix, which is not available in classical QTL mapping.

```
r2 = image(ctlres, against="phenotypes")
```

Reconstruct the network and write two sif files. One sif file contains the full network, the other 
holds the edge summary network.

```
CTLnetwork(ctlres)
```

We can use Cytoscape to visualize the created network (available from [www.cytoscape.org](http://www.cytoscape.org// "www.cytoscape.org") )

### Example Data and Formats

Example cross object in csvr format (link)[https://github.com/DannyArends/CTLmapping/tree/master/Rctl/tests "Gene expression data"]

Example matrixes in tab delim format (link)[https://github.com/DannyArends/CTLmapping/tree/master/D/test/data "Metabolite abundance data"]

Example matrixes in experimental qtab format (link)[https://github.com/DannyArends/CTLmapping/tree/master/D/test/data "Metabolite abundance qtab data"]

### Contributing and TODO

Want to contribute? Great!

See Rctl/inst/TODO.txt or submit a Github issue or pull request

Its also possible to just post comments on code / commits.
Or be a maintainer, and adopt a function

Danny Arends

### Disclaimer

Copyright (c) 2010-2013 GBIC: Danny Arends, Ritsert C Jansen, Pjotr Prins

