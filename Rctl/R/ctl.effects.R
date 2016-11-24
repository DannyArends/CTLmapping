#
# ctl.effects.R
#
# copyright (c) 2016-2020 - GBIC, Danny Arends, Pjotr Prins, Yang Li, and Ritsert C. Jansen
# last modified Nov, 2016
# first written Nov, 2016
# 
# Memory efficient multi-threaded wrapper around the correlation and chisquare code
#

ctlmarker <- function(marker, phenotypes, phe, n.markers, verbose) {
  require(ctl)
  mgt <- unique(marker)
  if(any(is.na(mgt))) mgt <- mgt[!is.na(mgt)]
  phenotypes[is.na(phenotypes)] <- -999

  nsamples <- c()
  rmatrix <- matrix(NA, length(mgt), ncol(phenotypes), dimnames=list(mgt, colnames(phenotypes)))

  for (genotype in mgt) {
    indx <- which(marker == genotype)
    phenotype <- phenotypes[indx, phe]
    if (length(phenotype) > 0) {
      res <- rep(0, ncol(phenotypes));
      result <- .C("R_correlation1toN", x = as.double(phenotype), y = as.double(unlist(phenotypes[indx,])), res = as.double(res), 
                                        as.integer(length(phenotype)), as.integer(ncol(phenotypes)), as.integer(1), # No double multithreading
                                        as.integer(verbose), PACKAGE="ctl")
      rmatrix[genotype, ] <- result$res
      nsamples <- c(nsamples, length(phenotype))
    }
  }
  res <- rep(0, ncol(phenotypes))
  result <- .C("R_chiSQN", nr = as.integer(length(mgt)), r = as.double(unlist(t(rmatrix))), res = as.double(res), 
                           phe = as.integer(-1), nsamples = as.integer(nsamples), nphe = as.integer(ncol(phenotypes)), PACKAGE="ctl")
  raw.p = pchisq(result$res, (ncol(phenotypes)-1), 0, FALSE)
  adj.p = unlist(lapply(n.markers * ncol(phenotypes) * raw.p, min, 1.0))
  return( list(cors = rmatrix, chisq = result$res, raw.p = raw.p, adj.p = adj.p, lod = -log10(adj.p)) )
}

ctleffects <- function(genotypes, phenotypes, phe = 1, nthreads = 1, verbose = TRUE) {
  cl <- makeCluster(getOption("cl.cores", nthreads))
  m <- split(genotypes, rep(1:ncol(genotypes), each = nrow(genotypes)))
  r <- parLapply(cl, m, ctlmarker, phenotypes, phe, ncol(genotypes), verbose)
  names(r) <- colnames(genotypes)
  stopCluster(cl)
  return(r)
}

test.ctlmarker <- function(){
  require(ctl)
  require(parallel)
  data(ath.metabolites)
  res1 <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, phenocol=1, parametric = TRUE)
  res <- ctleffects(ath.metab$genotypes, ath.metab$phenotypes)

  dcors <- t(matrix(unlist(lapply(res,"[", 2)), ncol(ath.metab$phenotypes), ncol(ath.metab$genotypes)))
  raw.p <- t(matrix(unlist(lapply(res,"[", 3)), ncol(ath.metab$phenotypes), ncol(ath.metab$genotypes)))
  adj.p <- t(matrix(unlist(lapply(res,"[", 4)), ncol(ath.metab$phenotypes), ncol(ath.metab$genotypes)))
  LOD <- t(matrix(unlist(lapply(res,"[", 5)), ncol(ath.metab$phenotypes), ncol(ath.metab$genotypes)))

  op <- par(mfrow = c(2,1))
  image(res1[[1]]$dcor)
  image(dcors)
}
