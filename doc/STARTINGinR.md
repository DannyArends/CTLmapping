## Starting in R

Load the library in the R interface by the following command (in R):

```
library(ctl)                            # Load the library
?ctl                                    # Show the help
```

### Examples

If you have data prepared for (R/qtl)[http://www.rqtl.org/], you can use the 
CTLscan.cross() function to directly scan your R/qtl cross object.

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

We can use Cytoscape to visualize the created network (available from [www.cytoscape.org](http://www.cytoscape.org/ "www.cytoscape.org") )

### Example Data and Formats

Example cross object in (R/qtl csvr format) in the /Rctl/tests/ directory

Example matrixes in (tab delim format) in the /D/test/data/ directory

Example matrixes in (experimental qtab format) in the /D/test/data/ directory

