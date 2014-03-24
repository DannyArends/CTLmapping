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
    if(nrow(y) != length(x)) stop("incompatible dimensions")
    res <- rep(0, ncol(y));
    result <- .C("R_correlation1toN", x = as.double(x), y = as.double(y), res = as.double(res), as.integer(length(x)), as.integer(ncol(y)), as.integer(verbose), PACKAGE="ctl")
  }else{
    if(length(y) != length(x)) stop("incompatible dimensions")
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

chiSQtoP <- function(cv, dof){
  unlist(lapply(cv,function(x){
    res <- 0;
    result <- .C("R_chiSQtoP", Cv = as.double(x), Dof = as.integer(dof), res = as.double(res), PACKAGE="ctl")
    return(result$res)
  }))
}

test.correlation <- function(){
  library(ctl)
  correlation(1:10, 1:10)
  correlation(10:1, 1:10)
  v <- runif(10)
  against <- matrix(runif(50),10,5)
  correlation(v, against)     # This
  cor(v,against)              # Should match
  correlation(v, t(against))  # Here we check the error condition
  cor(v,t(against))

}

test.chiSQN <- function(){
  library(ctl)
  ngeno <- 2
  correlations <- matrix(runif(20), 10, ngeno)  # 10 markers, with 2 x correlation per genotype
  nind <- c(20, 20)                             # 20 individuals in each group
  cvs <- chiSQN(correlations, nind)             # Calculate the critical ChiSquare values
  chiSQtoP(cvs, ngeno-1)                        # P values associated with the observed differences in correlation
}

# end of ctl.correlation.R


