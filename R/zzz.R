#
# qcl.significance.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified feb, 2011
# first written nov, 2010
# 
# .First.lib is run when the package is loaded with library(qtl)
#

.d_supported <- TRUE

.First.lib <- function(lib, pkg){
 library.dynam("QCL", pkg, lib)
 tryCatch(
   library.dynam("Dcode", pkg, lib),
   error = function(e){
     cat("No D code, falling back to R\n")
     .d_supported <<- FALSE
   })
}
# end of zzz.R
