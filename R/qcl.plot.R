#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Oct, 2011
# first written Nov, 2010
# 
# Plotting routines for QCL analysis
#

print.QCLscan <- function(x, ...){
  cat("QCLscan summary\n\n")
  cat("- Number of phenotypes scanned",length(x),"/",dim(x[[1]])[1],"\n")
  cat("- Number of markers",dim(x[[1]])[2],"\n")
  unlist(...)
}

plot.QCLscan <- function(x, cross, ...){
  npheno <- length(x)
  totpheno <- dim(x[[1]])[1]
  totmarkers <- dim(x[[1]])[2]
}
