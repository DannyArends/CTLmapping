#
# zzz.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified feb, 2011
# first written nov, 2010
# 
# .First.lib is run when the package is loaded with library(ctl)
#

#Global package variables
.CtlEnv <- new.env()

get_CtlEnv <- function(){.CtlEnv}
has_snow <- function(){ get(".has_snow", envir = .CtlEnv) }
has_rqtl <- function(){ get(".has_rqtl", envir = .CtlEnv) }

.has_snow_warnmsg    <- "- Package: SNOW not found, multi CPU support is not available/disabled\n"
.has_rqtl_warnmsg    <- "- Package: qtl not found, R/qtl functionality not available/disabled\n"

#Package loading
.onAttach <- function(lib, pkg){
  packageStartupMessage("- Loading package ctl\n", appendLF = FALSE)
  .has_snow <- FALSE
  .has_rqtl <- FALSE
  .has_snow  <- ("snow" %in% installed.packages(lib.loc=lib)[,1])
  .has_rqtl  <- ("qtl" %in% installed.packages(lib.loc=lib)[,1])
  assign(".has_snow",  .has_snow,  envir = .CtlEnv)
  assign(".has_rqtl",  .has_rqtl,  envir = .CtlEnv)
  if(!get(".has_snow", envir = .CtlEnv)){
    packageStartupMessage(.has_snow_warnmsg, appendLF = FALSE)
  }
  if(!get(".has_rqtl", envir = .CtlEnv)){
    packageStartupMessage(.has_rqtl_warnmsg, appendLF = FALSE)
  }
}

# end of zzz.R
