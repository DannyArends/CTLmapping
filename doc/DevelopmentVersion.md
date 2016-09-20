### Installing a development version

#### Option 1) Installing the development version
You can install a development version from the github repository by using the 
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

#### Option 2) Manually install the R source code

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

