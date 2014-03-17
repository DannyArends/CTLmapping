#
# ctl.correlation.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Pjotr Prins, Yang Li, and Ritsert C. Jansen
# last modified Dec, 2012
# first written Dec, 2012
# 
# Wrappers around the correlation code
#

correlation <- function(x, y, verbose = FALSE){
  res <- c(0);
  if(is.matrix(y)){
    res <- rep(0, ncol(y));
    result <- .C("R_correlation1toN", x = as.double(x), y = as.double(y), res = as.double(res), as.integer(length(x)), as.integer(ncol(y)), as.integer(verbose), PACKAGE="ctl")
  }else{
    result <- .C("R_correlation", x = as.double(x), y = as.double(y), res = as.double(res), as.integer(length(x)), as.integer(verbose), PACKAGE="ctl")
  }
  return(result$res)
}


#void R_chiSQN(int* nr, double* r, double* res, int* phe, int* nsamples, int* nphe){

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

# end of ctl.correlation.R

