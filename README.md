## Correlated Trait Locus (CTL) mapping

In this repository implementations of the novel Correlated Trait Locus (CTL) 
mapping algorithm. Provided as package for the statistical language R and in 
the D 2.0 programming language.

CTL mapping is a novel approach to detect genetic regulation of phenotypes in 
natural and experimental populations. It is a method which complements classical 
QTL analysis, providing additional insights overlooked by the classical QTL 
approach.

### Installing the R package
The quickest way to start mapping CTLs, is to install the package directly into 
R from Github using the devtools package:

```
# install.packages("devtools")
library(devtools)
install_github("CTLmapping", "DannyArends", subdir="Rctl")
```

Learn more about the [the R commands](https://github.com/DannyArends/CTLmapping/blob/master/Learn CTL/STARTINGinR.md)

### Download the software (R, C and D)

The second option to is to clone the package from Github, this will give you 
access to the R package but also provide the command line tools. First prepare 
your environment by download and 'moving' to the folder:

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping                                          # Goto the folder

### Manually install the R library

Prepare your environment by download and installing the R environment from 
[www.r-project.org](http://www.r-project.org/ "www.r-project.org"). Then 
download CTLmapping and install into R by using (from a terminal / command 
line):

    $ R CMD INSTALL Rctl                                     # Install the package

or use the 'installR' makefile target:

    $ make installR                                          # Install into R

Plans are to put the package on CRAN, but this has not happend yet. A Quick online 
introduction is [available](https://github.com/DannyArends/CTLmapping/blob/master/Learn CTL/STARTINGinR.md) 
but help files are also easily available for almost all function in R using:

```
library(ctl)                            # Load the library
?ctl                                    # Show the general help for ctl
```

### Compile the standalone executable

Optimized versions of the software are also available for more high throughput data.
Here we explain how to build a standalone executable version using C

#### (C version)

Just run 'make' from a terminal / command line:

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

    $ make versionD                                          # Compile the executable

Optionally you can use a provided binary by downloading the approriate one for your 
operating system.

### Algorithm

Learn more about the [algorithm](https://github.com/DannyArends/CTLmapping/blob/master/Learn CTL/ALGORITHM.md)

### Contributing and TODO

Want to contribute? Great!

See Rctl/inst/TODO.txt or submit a Github issue or pull request

Its also possible to just post comments on code / commits.
Or be a maintainer, and adopt a function

Danny Arends

### Disclaimer

Copyright (c) 2010-2013 GBIC: Danny Arends, Ritsert C Jansen, Pjotr Prins

