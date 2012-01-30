#
# qcl.significance.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified feb, 2011
# first written nov, 2010
# 
# .First.lib is run when the package is loaded with library(qtl)
#

#Global package variables
.QclEnv <- new.env()

get_QclEnv <- function(){.QclEnv}
has_snow <- function(){ get(".has_snow", envir = .QclEnv) }
has_rqtl <- function(){ get(".has_rqtl", envir = .QclEnv) }

.has_snow_warnmsg    <- "- Package: SNOW not found, multi CPU support is not available/disabled\n"
.has_rqtl_warnmsg    <- "- Package: qtl not found, R/qtl functionality not available/disabled\n"

#Package loading
.onAttach <- function(lib, pkg){
  packageStartupMessage("- Loading package qcl\n", appendLF = FALSE)
  .has_snow <- FALSE
  .has_rqtl <- FALSE
  .has_snow  <- ("snow" %in% installed.packages(lib.loc=lib)[,1])
  .has_rqtl  <- ("qtl" %in% installed.packages(lib.loc=lib)[,1])
  assign(".has_snow",  .has_snow,  envir = .QclEnv)
  assign(".has_rqtl",  .has_rqtl,  envir = .QclEnv)
  if(!get(".has_snow", envir = .QclEnv)){
    packageStartupMessage(.has_snow_warnmsg, appendLF = FALSE)
  }
  if(!get(".has_rqtl", envir = .QclEnv)){
    packageStartupMessage(.has_rqtl_warnmsg, appendLF = FALSE)
  }
}

# end of zzz.R
