#
# qcl.significance.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified feb, 2011
# first written nov, 2010
# 
# .First.lib is run when the package is loaded with library(qtl)
#
.First.lib <- function(lib, pkg){
 library.dynam("QCL", pkg, lib)
 library.dynam("D_QCL", pkg, lib)
}
# end of zzz.R
