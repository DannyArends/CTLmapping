mapCTL: Correlated Trait Locus (CTL) mapping in D
=================================================

Download the package
--------------------
Prepare your environment by download and 'moving' to the folder:

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping                                          # Goto the folder

Compile the D 2.0 version
-------------------------
Prepare your environment by download and installing the DMD 2.0 compiler from 
[www.d-programming-language.org](http://www.d-programming-language.org 
"www.d-programming-language.org"). Run 'make' from a terminal / command line:

    $ make versionD                                          # Compile the executable

Additional tools
----------------
Additionally we provide support for the QTAB file format available with 
qtlHD. To use QTAB files as input clone the qtlHD directory in the same 
folder as the CTLmapping repository and use rake to compile, the QTAB 
library is build and wrapped automatically:

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping                                          # Goto the folder
    $ make DqtlHD                                            # Compile the executable
    $
    $ ./mapctl -f=qtab -p=test/data/multi_phenotypes.qtab -g=test/data/multi

Commandline options
-------------------

    -[-h]elp        - Show the help file
    -[-v]erbose     - Verbose mode
    -[-o]verwrite   - Overwrite previous output files
    -[-n]perms      - Number of permutations
    -[-p]henotypes  - File containing phenotypes
    -[-g]enotypes   - File containing genotypes
    -[-f]ormat      - File format

Disclaimer
----------
Copyright (c) 2010 Danny Arends
