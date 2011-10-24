#
# QCL.R
#
# copyright (c) 2011 Danny Arends and Ritsert C. Jansen
# last modified Okt, 2011
# first written nov, 2010
# 
# D functions to do correlations

correlation <- function(x,y){
  if(is.matrix(x) && dim(x)[2] > 1){
    #warning("'x' is a matrix, ignoring 'y'")
    return(matrix_correlation(x));
  }else{
    if(missing(y))stop("'y' is missing");
    if(length(x) == length(y)){
      return(.C("correlation_v",
        x=as.numeric(x),
        y=as.numeric(y),
        as.integer(length(x)),
        cor=as.numeric(0),
        PACKAGE="qcl")$cor
      )
    }else{
      stop("The length of 'x' is not equal to the length of 'y'")
    }
  }
}

matrix_correlation <- function(x){
  res <- .C("correlation_m",
      nrx=as.integer(nrow(x)),
      ncx=as.integer(ncol(x)),
      x=as.numeric(x),
      cor=as.numeric(rep(1,ncol(x)*ncol(x))),
      PACKAGE="qcl"
    )
  mcor = matrix(res$cor,ncol(x),ncol(x))
  mcor
}
