#
# zzz.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Nov, 2010
# 
# .onLoad is run when the package is loaded with library(ctl)

.onLoad <- function(lib, pkg) library.dynam("ctl", pkg, lib)

# end of zzz.R
