## Correlated Trait Locus (CTL) mapping

In this repository implementations of the novel Correlated Trait Locus (CTL) 
mapping algorithm. Provided as package for the statistical language R and in 
the D 2.0 programming language.

CTL mapping is a novel approach to detect genetic regulation of phenotypes in 
natural and experimental populations. It is a method which complements classical 
QTL analysis, providing additional insights overlooked by the classical QTL 
approach.

#### Algorithm

Learn more about the [algorithm](https://github.com/DannyArends/CTLmapping/blob/master/learn%20CTL/ALGORITHM.md) 
behind CTL mapping

### Option 1) Installing the R package
The quickest and prefered way to start mapping CTLs, is to install the package directly into 
R from Github using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) 
package. On windows Rtools needs to be installed on the system. Download and 
install the version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/) that matches 
your current R version.

After installation, run the following commands in the R terminal to install the CTL mapping package.
Uncomment the first line (by removing the #), if you do not have devtools installed yet:

```
# install.packages("devtools")
library(devtools)
install_github("DannyArends/CTLmapping", subdir="Rctl")
```

After this learn more about the [the R commands](https://github.com/DannyArends/CTLmapping/blob/master/learn%20CTL/STARTINGinR.md) 
to start mapping CTLs on example data and how to prepare your own experimental data.

### Option 2) Manually install the R library

The second (more complex) option to is to clone the package from Github, this will give 
you access to the R source code and also allows you to build the standalone linux executables.
On windows Rtools needs to be installed on your system. Download and install the 
version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/) that matches your 
current R version.

#### Download the software (R, C and D)

Prepare your environment by download and installing the R environment from [www.r-project.org](http://www.r-project.org/ "www.r-project.org"). 
Get the CTL mapping source code by download and afterwards 'move' into the folder:

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping                                          # Goto the folder

The downloaded CTLmapping can be installed into R by using the following command (from a terminal / command line):

    $ R CMD INSTALL Rctl                                     # Install the package

or use the 'installR' makefile target:

    $ make installR                                          # Install into R

Plans are to put the package on CRAN, but this has not happend yet.

### Getting help in R
A quick online introduction is [available](https://github.com/DannyArends/CTLmapping/blob/master/learn%20CTL/STARTINGinR.md) 
but help files are also easily available for almost all function in R using:

```
library(ctl)                            # Load the library
?ctl                                    # Show the general help for ctl
?CTLscan                                # Show the help for the CTLscan function
```

### Compile the standalone executable
Note: makefiles to produce standalone executable are linux specific and will now work under windows

Optimized versions of the software are also available for more high throughput data. 
The follwoing two paragraphs explain how to build a standalone executable version 
using either the C front-end code or the D 2.0 front-end code.

#### (C version)

Clone the repository from github, move into the CTLmappinf folder and run 'make' from a terminal / command line:

    $ make versionC                                          # Compile the executable
    $ make static                                            # Compile the static library
    $ make shared                                            # Compile the shared library

C code can also be compiled also into static or dynamic libraries, when using the 
appropriate makefile commands and can then be used in your own software as external 
dependancy with minimal coupling to your own software.

#### (D 2.0 version)

Prepare your environment by download and installing the DMD 2.0 compiler from 
[www.d-programming-language.org](http://www.d-programming-language.org 
"www.d-programming-language.org"). Run 'make' from a terminal / command line:

    $ make versionD                                          # Compile the D 2.0 executable

In the future, you can use a provided binary by downloading the approriate one for your 
operating system.

### Issues

Issues can be raised through the github issue tracker.

### Tests

CTL mapping uses the build in R framework to test the package for global regressions and unit-testing of documented functions.
Tests can be executed from the commandline, by using the following command:

    $ R CMD check Rctl                                     # Run the unit-tests of the R package

### Contributing 

Want to contribute? Great! We're actively looking for someone to do the website 
[www.mapctl.org](http://www.mapctl.org "www.mapctl.org")

Contribute to CTL source code by forking the Github repository, and sending us pull requests.

see Rctl/inst/TODO.txt or submit a Github issue or pull request

Its also possible to just post comments on code / commits.
Or be a maintainer, and adopt a function

### Cite CTL mapping

When this software was helpful, please cite the JOSS paper: 
[TODO: http://joss.theoj.org/papers/ DOI] in your publication, when (re-)using the any 
code or data provided by this package, please also cite: [TODO: Github/Zenodo DOI] 
for future code and data references. Citations are also available for import in 
bibtex format [TODO: bibtex with DOIs].

### License

The CTL mapping source code is released under the GNU GENERAL PUBLIC LICENSE Version 3 (GPLv3). See LICENSE.txt.

This software was developed between 2012-2016 at the Groningen Bioinformatics Centre by Danny Arends, Yang Li, Pjotr Prins and Ritsert C. Jansen
