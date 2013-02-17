## mapCTL: Correlated Trait Locus (CTL) mapping in D

Prepare your environment by download and 'moving' to the folder:

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping/D                                        # Go to the D folder

### Compile the D 2.0 version

Prepare your environment by download and installing the DMD 2.0 compiler from 
[www.d-programming-language.org](http://www.d-programming-language.org 
"www.d-programming-language.org"). Run 'make' from a terminal / command line:

    $ make                                                   # Compile the executable

### Additional tools

Additionally we provide support for the QTAB file format available with 
qtlHD. To use QTAB files as input clone the qtlHD directory in the same 
folder as the CTLmapping repository and use rake to compile, the QTAB 
library is build and wrapped automatically:

    $ git clone git://github.com/DannyArends/qtlHD.git       # Download the qtlHD repository
    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping                                          # Goto the folder
    $ make DqtlHD                                            # Compile the executable
    $
    $ ./mapctl -f=qtab -p=test/data/multi_phenotypes.qtab -g=test/data/multi

### Commandline options

    -[-h]elp        - Show the help file
    -[-v]erbose     - Verbose mode

    -[-p]henotypes  - File containing phenotypes
    -[-g]enotypes   - File containing genotypes
    -[-f]ormat      - File format

    -[-o]ut         - Name of the output folder
    -[-r]edo        - Overwrite previous output files
    -[-n]perms      - Number of permutations

    --alpha         - Set Alpha (Default: 1)
    --beta          - Set Beta (Default: 1)

### Disclaimer

Copyright (c) 2010-2012 GBIC: Danny Arends, Ritsert C Jansen, Pjotr Prins

