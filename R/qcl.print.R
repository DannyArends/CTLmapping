#
# qcl.print.R
#
# copyright (c) 2010-2011 Danny Arends and Ritsert C. Jansen
# last modified Jan, 2012
# first written Oct, 2011
# 
# Print routines for QCL analysis
#

print.QCLscan <- function(x, ...){
  cat("QCLscan summary",attr(x$qcl,"name"),"\n\n")
  cat("- Number of background phenotypes",dim(x$qcl)[1],"\n")
  cat("- Number of markers",dim(x$qcl)[2],"\n")
  cat("- Number of permutations",length(unlist(x$p)),"\n")
  getPermuteThresholds(x,..., verbose = TRUE)
}

getPermuteThresholds <-function(x, ..., verbose = FALSE){
  if(any(class(x)=="QCLscan")){
    n <- dim(x$qcl)[2]
    x <- x$p
    if(!any(class(x)=="QCLpermute")) stop("No permutations found in the QCLscan object")
  }else{
    stop("No QCLscan object")
  }
  sorted <- sort(unlist(x))
  l <- length(sorted)
  values <- NULL
  valnames <- NULL
  for(x in c(.95,.99,.999)){
    if(1/((1-x)/n) < length(sorted)){
    v <- sorted[l*(1-(1-x)/n)]
    values <- c(values,v)
    valnames <- c(valnames,paste((1-x)*100,"%"))
    if(verbose) cat((1-x)*100,"%\t",v,"\n")
    }else{
    values <- c(values,NaN)
    valnames <- c(valnames,paste((1-x)*100,"%"))
    if(verbose) cat((1-x)*100,"%\t",NaN,"\n")
    }
  }
  names(values) <- valnames
  values
}

print.QCLpermute <- function(x, ...){
  getPermuteThresholds(x, ..., verbose=TRUE)
}

# end of qcl.print.R
