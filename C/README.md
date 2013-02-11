mapCTL: Correlated Trait Locus (CTL) mapping in C
=================================================

Download the package
--------------------
Prepare your environment by download and 'moving' to the folder:

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping                                          # Goto the folder

Compile the standalone
-------------------------
Prepare your environment by download and installing the DMD 2.0 compiler from 
[www.d-programming-language.org](http://www.d-programming-language.org 
"www.d-programming-language.org"). Run 'make' from a terminal / command line:

    $ make                                                   # Compile the executable
    $ make static                                            # Compile the static library
    $ make shared                                            # Compile the shared library

Commandline options
-------------------

    -h          - Show the help file

    -p          - File containing phenotypes
    -g          - File containing genotypes

    -n          - Number of permutations

    -a          - Set Alpha (Default: 1)
    -b          - Set Beta (Default: 1)

Disclaimer
----------
Copyright (c) 2010-2012 GBIC: Danny Arends, Ritsert C Jansen, Pjotr Prins

