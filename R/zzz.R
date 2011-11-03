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
assign(".has_d",     TRUE, envir = .QclEnv)
assign(".has_snow",  FALSE, envir = .QclEnv)
assign(".has_rcurl", FALSE, envir = .QclEnv)
assign(".has_rqtl",  FALSE, envir = .QclEnv)


.has_d_warnmsg       <- "- D not available, using standard R/C/C++ functionality\n"
.has_snow_warnmsg    <- "- Package: SNOW not found, multi CPU support is not available/disabled\n"
.has_rcurl_warnmsg   <- "- Package: RCURL not found, automatic KEGG pathway annotations not available/disabled\n"
.has_rqtl_warnmsg    <- "- Package: qtl not found, R/qtl functionality not available/disabled\n"

#Package loading
.onAttach <- function(lib, pkg){
  packageStartupMessage("- Loading package qcl\n")
  .has_d <- TRUE
  .has_snow <- FALSE
  .has_rcurl <- FALSE
  .has_rqtl <- FALSE
  library.dynam("qcl", pkg, lib)
  tryCatch(
    library.dynam("Dcode", pkg, lib),
    error = function(e){
     .has_d <<- FALSE
   })
  .has_snow  <<- ("snow" %in% installed.packages()[,1])
  .has_rcurl <<- ("RCurl" %in% installed.packages()[,1])
  .has_rqtl  <<- ("qtl" %in% installed.packages()[,1])
  assign(".has_d",     .has_d,     envir = .QclEnv)
  assign(".has_snow",  .has_snow,  envir = .QclEnv)
  assign(".has_rcurl", .has_rcurl, envir = .QclEnv)
  assign(".has_rqtl",  .has_rqtl,  envir = .QclEnv)
  if(!get(".has_d", envir = .QclEnv)) packageStartupMessage(.has_d_warnmsg)
  if(!get(".has_snow", envir = .QclEnv)) packageStartupMessage(.has_snow_warnmsg)
  if(!get(".has_rcurl", envir = .QclEnv)) packageStartupMessage(.has_rcurl_warnmsg)
  if(!get(".has_rqtl", envir = .QclEnv)) packageStartupMessage(.has_rqtl_warnmsg)
}

# end of zzz.R
