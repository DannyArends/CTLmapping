## Correlated Trait Locus (CTL) mapping

This repository contains the source code for the Correlated Trait Locus (CTL) mapping algorithm. 
CTL mapping is a novel approach to detect genetic regulation of phenotypes in natural and 
experimental populations. It is a method which complements classical QTL analysis, providing 
additional insights overlooked by the classical QTL approach.

#### Algorithm

Learn more about the [algorithm](https://github.com/DannyArends/CTLmapping/blob/master/learn%20CTL/ALGORITHM.md) 
behind CTL mapping

### Install

#### Option 1) Installing the R package from CRAN
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

#### Option 2) Installing the development version
Sometimes it is needed to install a develpoment version, since CRAN 
might take a while to update after a bug is fixed. You can install 
a development version from the github repository by using the 
[devtools](https://cran.r-project.org/web/packages/devtools/index.html) 
package. On windows Rtools needs to be installed on the system. So first download and 
install the version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/) that matches 
your current R version.

After installation of Rtools, run the following commands in the R 
terminal to download and install the CTL mapping package. Uncomment 
the first line (by removing the #), if you do not have devtools installed yet:

```R
# install.packages("devtools")
library(devtools)                                         # Load the devtools package
install_github("DannyArends/CTLmapping", subdir="Rctl")   # Install the package from Github
```

After this learn more about the [the R commands](learn%20CTL/STARTINGinR.md) 
to start mapping CTLs on example data and how to prepare your own experimental data.

#### Option 3) Manually install the R source code

The second (more complex) option to is to clone the package from Github, this will give 
you access to the R source code and also allows you to build the [standalone linux executables](learn%20CTL/COMPILE.md) .
On windows Rtools needs to be installed on your system. Download and install the 
version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/) that matches your 
current R version.

Prepare your environment by download and installing the R environment from [www.r-project.org](http://www.r-project.org/ "www.r-project.org"). 
Get the CTL mapping source code by download and afterwards 'move' into the folder:

    git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    cd CTLmapping                                          # Goto the folder

The downloaded CTLmapping can be installed into R by using the following command (from a terminal / command line):

    R CMD INSTALL Rctl                                     # Install the package

or use the 'installR' makefile target:

    make installR                                          # Install into R

Plans are to put the package on CRAN, but this has not happend yet.

### Test

CTL mapping uses the build in R framework to test the package for global regressions and unit-testing of documented functions.
Tests can be executed from the commandline, by using the following command:

    R CMD check Rctl                                      # Run the unit-tests of the R package

### Documentation
A short online introduction is [available](learn%20CTL/STARTINGinR.md) 
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

When this software was helpful, please cite the JOSS paper: 
[TODO: http://joss.theoj.org/papers/ DOI] in your publication, when (re-)using the any 
code or data provided by this package, please also cite: 

[![DOI](https://zenodo.org/badge/2619708.svg)](https://zenodo.org/badge/latestdoi/2619708)
 
Citations are also available for import in bibtex format [TODO: bibtex with DOIs].

### License

The CTL mapping source code is released under the GNU GENERAL PUBLIC LICENSE Version 3 (GPLv3). See [LICENSE.txt](LICENSE.txt).

This software was developed between 2012-2016 at the Groningen Bioinformatics Centre by Danny Arends, Yang Li, Pjotr Prins and Ritsert C. Jansen

### Contact

Code managed by Dr. Danny Arends and the Groningen Bioinformatics Centre, Groningen, NLD. 