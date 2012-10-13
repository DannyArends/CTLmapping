#
# zzz.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified feb, 2011
# first written nov, 2010
# 
# .onLoad is run when the package is loaded with library(ctl)
#

has_snow <- function() any(grep("snow",rownames(installed.packages())))
has_rqtl <- function() any(grep("qtl",rownames(installed.packages())))


.onLoad <- function(lib, pkg) library.dynam("ctl", pkg, lib)

# end of zzz.R
