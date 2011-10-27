#
# correlations.R
#
# copyright (c) 2011 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Oct, 2011
# 
# R->D functions to do correlations
#

correlation <- function(x,y){
  if(!.d_supported){
   return(cor(x,y))
  }else{
    if(is.matrix(x) && dim(x)[2] > 1){
      return(matrix_correlation(x))
    }else{
      if(missing(y))stop("'y' is missing");
      if(length(x) == length(y)){
        return(.C("correlation_v",
          x=as.numeric(x),
          y=as.numeric(y),
          as.integer(length(x)),
          cor=as.numeric(0),
          PACKAGE="Dcode")$cor
        )
      }else{
        stop("The length of 'x' is not equal to the length of 'y'")
      }
    }
  }
}

matrix_correlation <- function(x){
  if(!.d_supported){
    return(cor(x))
  }else{
    res <- .C("correlation_m",
      nrx=as.integer(nrow(x)),
      ncx=as.integer(ncol(x)),
      x=as.numeric(x),
      cor=as.numeric(rep(1,ncol(x)*ncol(x))),
      PACKAGE="Dcode"
    )
    mcor = matrix(res$cor,ncol(x),ncol(x))
    return(mcor)
  }
}
