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
.has_d               <- TRUE
.has_snow            <- FALSE
.has_rcurl           <- FALSE
.has_rqtl            <- FALSE

.has_d_warnmsg       <- "- D not available, using standard R/C/C++ functionality\n"
.has_snow_warnmsg    <- "- Package: SNOW not found, multi CPU support is not available/disabled\n"
.has_rcurl_warnmsg   <- "- Package: RCURL not found, automatic KEGG pathway annotations not available/disabled\n"
.has_rqtl_warnmsg    <- "- Package: qtl not found, R/qtl functionality not available/disabled\n"

#Package loading
.First.lib <- function(lib, pkg){
 cat("- Loading package qcl\n")
 library.dynam("qcl", pkg, lib)
 tryCatch(
   library.dynam("Dcode", pkg, lib),
   error = function(e){
     .has_d <<- FALSE
   })
  .has_snow  <<- ("snow" %in% installed.packages()[,1])
  .has_rcurl <<- ("RCurl" %in% installed.packages()[,1])
  .has_rqtl  <<- ("qtl" %in% installed.packages()[,1])
  
  if(!.has_d)     cat(.has_d_warnmsg)
  if(!.has_snow)  cat(.has_snow_warnmsg)
  if(!.has_rcurl) cat(.has_rcurl_warnmsg)
  if(!.has_rqtl)  cat(.has_rqtl_warnmsg)
}

# end of zzz.R
