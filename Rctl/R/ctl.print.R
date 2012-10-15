#
# ctl.print.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Oct, 2011
# 
# Print routines for CTL analysis
#

print.CTLobject <- function(x, ...){
  cat("CTLobject summary\n\n")
  cat("- Number of scanned phenotypes:",length(x),"\n")
}

print.CTLscan <- function(x, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  cat("CTLscan summary",attr(x$ctl,"name"),"\n\n")
  cat("- Number of background phenotypes",dim(x$ctl)[1],"\n")
  cat("- Number of markers",dim(x$ctl)[2],"\n")
  cat("- Number of permutations",length(unlist(x$p)),"\n")
  invisible(CTLsignificant(x))
}

getPermuteThresholds <-function(x, significance = c(.05,.01,.001), ..., verbose = FALSE){
  sorted <- sort(unlist(x))
  l <- length(sorted)
  values <- NULL
  valnames <- NULL
  for(x in significance){
    if(1/(x) < length(sorted)){
      v <- sorted[l*(1-(x))]
      values <- c(values,v)
      valnames <- c(valnames,paste(x*100,"%"))
      if(verbose) cat(x*100,"%\t",v,"\n")
    }else{
      values <- c(values,NaN)
      valnames <- c(valnames,paste(x*100,"%"))
      if(verbose) cat(x*100,"%\t",NaN,"\n")
    }
  }
  names(values) <- valnames
  values
}

print.CTLpermute <- function(x, ...){
  getPermuteThresholds(x, ..., verbose=TRUE)
}

# end of ctl.print.R
