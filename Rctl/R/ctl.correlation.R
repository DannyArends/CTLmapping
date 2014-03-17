#
# ctl.correlation.R
#
# copyright (c) 2010-2014 - GBIC, Danny Arends, Pjotr Prins, Yang Li, and Ritsert C. Jansen
# last modified Mar, 2014
# first written Mar, 2014
# 
# Wrappers around the correlation and chisquare code
#

correlation <- function(x, y, verbose = FALSE){
  if(is.matrix(y)){           # Call the optimized loop unrolled version
    res <- rep(0, ncol(y));
    result <- .C("R_correlation1toN", x = as.double(x), y = as.double(y), res = as.double(res), as.integer(length(x)), as.integer(ncol(y)), as.integer(verbose), PACKAGE="ctl")
  }else{
    res <- c(0);
    result <- .C("R_correlation", x = as.double(x), y = as.double(y), res = as.double(res), as.integer(length(x)), as.integer(verbose), PACKAGE="ctl")
  }
  return(result$res)
}

chiSQN <- function(correlations, nsamples){
  res <- rep(0, nrow(correlations));
  result <- .C("R_chiSQN", nr = as.integer(ncol(correlations)), as.double(correlations), res = as.double(res), as.integer(-1), as.integer(nsamples), nrow(correlations), PACKAGE="ctl")
  return(result$res)
}

test.correlation <- function(){
  library(ctl)
  correlation(1:10, 1:10)
  correlation(10:1, 1:10)
  correlation(runif(10), matrix(runif(50),10,5))
}

test.chiSQN <- function(){
  library(ctl)
  correlations <- matrix(runif(6), 3, 2)
  nind <- c(20,20)
  chiSQN(correlations, nind)
}

# end of ctl.correlation.R

