## mapCTL: Correlated Trait Locus (CTL) mapping in C

Prepare your environment by download and 'moving' to the folder:

    $ git clone git://github.com/DannyArends/CTLmapping.git  # Download the repository
    $ cd CTLmapping/C                                        # Go to the C folder

### Compilation

Prepare your environment by download and installing a C compiler. Run 'make' from a 
terminal / command line to build either the executable, shared or static library:

    $ make R_HOME="path/to/R"               # Compile the executable
    $ make R_HOME="path/to/R" static        # Compile the static library
    $ make R_HOME="path/to/R" shared        # Compile the shared library

### Commandline options

The executable supports the following command line options:

    -p<FILE>   Input file with phenotype data (Default: phenotypes.csv)
    -g<FILE>   Input file with genotype data (Default: genotypes.csv)
    -o<FILE>   Name of the output file (Default: summary.txt)
    -t<N>      Significance threshold (0..1) (Default: 0.01)
    -f         Save all the results to file
    -d         When specified permutation are performed
    -n<N>      # of permutations (Default: 100)
    -h         Shows this help

### Testing the executable

To test the executable run a special make target:

    $ make test                                              # Compiles a test executable
    $ ./mapctl

### Disclaimer

(c) 2012-2020 GBIC - RUG & HU-Berlin, written by Danny Arends, Ritsert C Jansen, Pjotr Prins
