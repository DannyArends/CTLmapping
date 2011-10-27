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
.d_supported  <- TRUE
.d_warningmsg <- "- D not available, using standard R/C/C++ functionality\n"

#Package loading
.First.lib <- function(lib, pkg){
 cat("- Loading package QCL\n")
 library.dynam("QCL", pkg, lib)
 tryCatch(
   library.dynam("Dcode", pkg, lib),
   error = function(e){
     cat(.d_warningmsg)
     .d_supported <<- FALSE
   })
}

# end of zzz.R
