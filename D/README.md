mapCTL: Correlated Trait Locus (CTL) mapping in D
=================================================

Installation
------------
Prepare your environment by following these steps:

- Download and Install the DMD compiler from [www.digitalmars.com/d/download.html](http://www.digitalmars.com/d/download.html "www.digitalmars.com/d/download.html")
- Install Ruby 1.9.1 (or higher) from (http://www.ruby-lang.org/)
- Install Rake from (http://rake.rubyforge.org/)
- Clone the repository

```
    git clone https://DannyArends@github.com/DannyArends/CTLmapping.git
    cd CTLmapping/D
    compile
    ./mapctl
```

Additional tools
----------------
Additionally we provide support for the QTAB file format available with 
qtlHD. To use QTAB files as input clone the qtlHD directory in the same 
folder as the CTLmapping repository and use rake to compile, the QTAB 
library is build and wrapped automatically:

```
    git clone https://DannyArends@github.com/DannyArends/CTLmapping.git
    git clone https://DannyArends@github.com/DannyArends/qtlHD.git
    cd CTLmapping/D
    rake
    ./mapctl -f=qtab -p=test/data/multi_phenotypes.qtab -g=test/data/multi
```

Single marker QTL mapping is also available by supplying the 'qtl' build 
target to rake:

```
    rake qtl
```

Commandline options
-------------------
```
    -[-h]elp        - Show the help file
    -[-v]erbose     - Verbose mode
    -[-o]verwrite   - Overwrite previous output files
    -[-n]perms      - Number of permutations
    -[-p]henotypes  - File containing phenotypes
    -[-g]enotypes   - File containing genotypes
    -[-f]ormat      - File format
```

Disclaimer
----------
Copyright (c) 2010 Danny Arends
