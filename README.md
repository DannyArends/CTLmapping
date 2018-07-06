## Correlated Trait Locus (CTL) mapping
[![Build Status](https://travis-ci.org/DannyArends/CTLmapping.svg)](https://travis-ci.org/DannyArends/CTLmapping)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ctl)](http://cran.r-project.org/package=ctl)
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.00087/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00087)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.154465.svg)](https://doi.org/10.5281/zenodo.154465)

This repository contains the source code for the Correlated Trait Locus (CTL) mapping algorithm. 
CTL mapping is a novel approach to detect genetic regulation of phenotypes in natural and 
experimental populations. It is a method which complements classical QTL analysis, providing 
additional insights overlooked by the classical QTL approach.

#### Algorithm

Learn more about the [algorithm](doc/ALGORITHM.md) behind CTL mapping. 
Alternatively I gave a workshop during the complex traits consortium 
meeting in Glasgow to demonstrate CTL mapping on an example data set. 
The presentation and code can be found on my [website](https://www.dannyarends.nl/ctc/)

### Installing the R package from CRAN
The quickest and prefered way to start mapping CTLs, is to 
install the package directly from CRAN, start R and issue 
the following command to install the package:

```R
install.packages("ctl")
```

After installation, load the package by:

```R
library("ctl")
```

Next thing is to [map your first CTL](doc/STARTINGinR.md) 
on either example data or your own data. For this see: [doc/STARTINGinR.md](doc/STARTINGinR.md) 

Sometimes it is needed to install a develpoment version, since CRAN 
might take a while to update after a bug is fixed. To learn how to 
install a development version see: [doc/DEVELOPMENT.md](doc/DEVELOPMENT.md) 

### Test

CTL mapping uses the build in R framework to test the package 
for global regressions and unit-testing of documented functions.
Tests can be executed from the commandline, by using the following command:

    R CMD check Rctl                                      # Run the unit-tests of the R package

### Documentation
A short online introduction is [available](doc/STARTINGinR.md) 
and help files with examples are also available for almost all functions in R using:

```R
library(ctl)                            # Load the library
?ctl                                    # Show the general help for ctl
?CTLscan                                # Show the help for the CTLscan function
```

### Issues

Issues can be raised through the github issue tracker.

### Contributing 

Want to contribute? Great! We're actively looking for someone to do the website 
[www.mapctl.org](http://www.mapctl.org "www.mapctl.org")

Contribute to CTL source code by forking the Github repository, and sending us pull requests.

For a list of active developments tasks, see [Rctl/inst/TODO.txt](Rctl/inst/TODO.txt) 

Its also possible to just post comments on code / commits.

Or be a maintainer, and adopt a function

### Cite

When this software was helpful, please cite the JOSS paper in your publication: 
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.00087/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00087)

When (re-)using the any code or data provided by this package, please also cite: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.154465.svg)](https://doi.org/10.5281/zenodo.154465)
 
Citations are also available for import in bibtex format [TODO: bibtex with DOIs].

### License

The CTL mapping source code is released under the GNU GENERAL PUBLIC LICENSE Version 3 (GPLv3). See [LICENSE.txt](LICENSE.txt).

This software was developed between 2012-2016 at the Groningen Bioinformatics Centre by Danny Arends, Yang Li, Pjotr Prins and Ritsert C. Jansen

### Contact

Code managed by Dr. Danny Arends and the Groningen Bioinformatics Centre, Groningen, NLD. 