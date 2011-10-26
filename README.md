QCL mapping
===========
QCL mapping is a novel approach to detect genetic regulation of phenotypes in recombinant 
inbred line populations (RIL). It is a method complementair to QTL analysis, and provides
additional insights, overlooked by the classical QTL approach  
 

Dependencies
------------
R software environment from [www.r-project.org](http://www.r-project.org/ "www.r-project.org")

Installation
------------
Prepare your environment by following these steps:

- Download and Install the R environment

Then install into R by using (from a terminal / commanline):

    $ git clone git://github.com/DannyArends/QCLmapping.git  # Download the repository
    $ R CMD INSTALL QCLmapping                               # Install the package

Optionally you can install the pre-build packages by downloading the appropriate 
package for your operating system. 

Starting
--------
Load the library in the R interface by the following command (in R):
    
    $ > library(qcl)                                # Load the library
    $ > ?qcl                                        # Show the help

example:
    
    $ > ?QCL.scan                                   # Show the help
    $ > data(multitrait)
    $ > res <- QCL.scan(multitrait)
    $ > res

QCL TODO
--------------------
See inst/TODO.txt

Contributing
------------

Want to contribute? Great!
Just post comments on code / commits.

Danny Arends

Disclaimer
----------
Copyright (c) 2010 Danny Arends
