#
# qcl.print.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Oct, 2011
# first written Oct, 2011
# 
# Plotting routines for QCL analysis
#

print.QCLscan <- function(x, ...){
  cat("QCLscan summary\n\n")
  cat("- Number of phenotypes scanned",length(x),"/",dim(x[[1]])[1],"\n")
  cat("- Number of markers",dim(x[[1]])[2],"\n")
}

print.QCLpermute <- function(x, ...){
  maximums <- lapply(x,function(x){apply(x,2,max)})
  sorted <- sort(unlist(maximums))
  l <- length(sorted)
  values <- NULL
  valnames <- NULL
  for(x in c(.95,.99,.999)){
    values <- c(values,sorted[l*x])
    valnames <- c(valnames,paste((1-x)*100,"%"))
    cat((1-x)*100,"%\t",sorted[l*x],"\n")
  }
  names(values) <- valnames
  invisible(values)
}
