#
# ctl.correlation.R
#
# copyright (c) 2010-2014 - GBIC, Danny Arends, Pjotr Prins, Yang Li, and Ritsert C. Jansen
# last modified Mar, 2014
# first written Mar, 2014
# 
# Wrappers around the correlation and chisquare code
#

correlation <- function(x, y, nthreads = 1, verbose = FALSE){
  if(is.matrix(y)){           # Call the optimized loop unrolled version
    if(nrow(y) != length(x)) stop(paste0("incompatible dimensions", length(x), nrow(y)))
    res <- rep(0, ncol(y));
    result <- .C("R_correlation1toN", x = as.double(x), y = as.double(y), res = as.double(res), 
                                      as.integer(length(x)), as.integer(ncol(y)), as.integer(nthreads),
                                      as.integer(verbose), PACKAGE="ctl")
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

getCorrelations <- function(genotypes, phenotypes, phenocol = 1, marker.col = 1, parametric = FALSE, verbose = TRUE){
  phenotypes <- apply(phenotypes, 2, as.numeric)
  marker <- genotypes[,marker.col]
  result <- NULL; nums <- NULL
  if(!parametric){
    st <- proc.time()
    phenotypes <- apply(phenotypes, 2, rank, na.last="keep") # Parametric vs Non-Parametric testing
    if(verbose) cat("Data ranking finished after:",(proc.time()-st)[3],"seconds\n")
  }
  for(x in names(table(marker))){
    idx <- which(marker == x)
    phenotypes[is.na(phenotypes)] <- -999
    pheX <- phenotypes[idx, phenocol]
    pheM <- phenotypes[idx, ]
    result$correlations <- cbind(result$correlations, correlation(pheX, pheM))
    result$samplesize <- c(result$samplesize, length(idx))
  }
  rownames(result$correlations) <- colnames(phenotypes)
  colnames(result$correlations) <- names(table(marker))
  return(result)
}

getCorrelations.cross <- function(cross, phenocol = 1, marker.col = 1){
  marker <- pull.geno(cross)[,marker.col]
  phenotypes <- pull.pheno(cross)
  result <- NULL; nums <- NULL
  for(x in names(table(marker))){
    idx <- which(marker == x)
    result$correlations <- cbind(result$correlations, correlation(phenotypes[idx, phenocol], phenotypes[idx, ]))
    result$samplesize <- c(result$samplesize, length(idx))
  }
  rownames(result$correlations) <- colnames(phenotypes)
  colnames(result$correlations) <- names(table(marker))
  return(result)
}

test.getCorrelation.cross <- function(){
  data(hyper, envir = environment())
  hyper$pheno <- matrix(runif(250*10), 250, 10)  # Create 10 phenotypes as a matrix
  colnames(hyper$pheno) <- paste0("Pheno", 1:10)
  crs <- getCorrelations.cross(hyper)
  cvs <- chiSQN(crs$correlations, crs$samplesize)
  pvs <- chiSQtoP(cvs, (length(crs$samplesize)-1))
}

test.correlation <- function(){
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
  ngeno <- 2
  correlations <- matrix(runif(20), 10, ngeno)  # 10 markers, with 2 x correlation per genotype
  nind <- c(20, 20)                             # 20 individuals in each group
  cvs <- chiSQN(correlations, nind)             # Calculate the critical ChiSquare values
  chiSQtoP(cvs, ngeno-1)                        # P values associated with the observed differences in correlation
  plot(chiSQtoP(seq(0,20,0.01), 1)  , lwd=2, t='l')                   # Plot the Pvalues for scores 0 till 30 using 1 degree of freedom
  points(chiSQtoP(seq(0,20,0.01), 2), lwd=2, t='l', col='yellow')     # Plot the Pvalues for scores 0 till 30 using 2 degree of freedom
  points(chiSQtoP(seq(0,20,0.01), 3), lwd=2, t='l', col='red')        # Plot the Pvalues for scores 0 till 30 using 3 degree of freedom
}

# end of ctl.correlation.R

