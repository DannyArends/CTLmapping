#
# ctl.slope.R
#
# copyright (c) 2010-2013 - GBIC, Danny Arends and Ritsert C. Jansen
# last modified Apr, 2013
# first written Apr, 2013
# 
# Slope analysis for CTL analysis
#

stdSlopeEstimate <- function(t1, t2, geno){
  betas <- NULL
  for(g in munique(geno)){
    ind    <- which(geno == g)
    t1mean <- mean(t1[ind], na.rm=TRUE)
    t2mean <- mean(t2[ind], na.rm=TRUE)
    nom <- 0; denom <- 0
    for(x in ind){
      if(!is.na(t1[x]) && !is.na(t2[x])){
        nom   <- nom + (t1[x] - t1mean) * (t2[x] - t2mean)
        denom <- denom + ((t1[x] - t1mean)^2)
      }
    }
    betas <- c(betas, nom/denom)
  }
  return(betas)
}

stdSlopeIntercept <- function(t1, t2, geno, betas){
  intercepts <- NULL
  group <- 1
  for(g in munique(geno)){
    ind    <- which(geno == g)
    t1mean <- mean(t1[ind], na.rm=TRUE)
    t2mean <- mean(t2[ind], na.rm=TRUE)
    intercepts <- c(intercepts, t2mean - betas[group] * t1mean)
    group <- group+1
  }
  return(intercepts)
}

stdSlopeError <- function(t1, t2, geno, betas){
  inters <- stdSlopeIntercept(t1, t2, geno, betas)
  errors <- NULL
  group <- 1
  for(g in munique(geno)){
    ind    <- which(geno == g)
    Nj <- length(ind)
    error  <- 0
    for(x in ind){
      if(!is.na(t1[x]) && !is.na(t2[x])){
        yhat  <- inters[group] + (betas[group] * t1[x])
        error <- error + (t2[x] - yhat)^2
      }
    }
    error <- error / (Nj-2)
    errors <- c(errors, error)
    group <- group + 1
  }
  return(errors)
}

stdSlopeBeta <- function(t1, t2, geno, betas){
  nom <- 0; denom <- 0; group <- 1

  for(g in munique(geno)){
    ind <- which(geno == g)
    t1mean <- mean(t1[ind], na.rm=TRUE)
    value <- 0
    for(x in ind){ 
      if(!is.na(t1[x])){ value <- value + (t1[x] - t1mean)^2 }
    }
    nom <- nom + (betas[group] * value)
    denom <- denom + value
    group <- group+1
  }
  return(nom/denom)
}

stdSlopeTest <- function(t1, t2, geno, verbose = FALSE){
  betas   <- stdSlopeEstimate(t1, t2, geno)       # Slope estimates per genotype
  errors  <- stdSlopeError(t1, t2, geno, betas)   # Squared error
  stdbeta <- stdSlopeBeta(t1, t2, geno, betas)    # Slope on all data
  if(verbose){ cat("Betas", betas, "\n"); cat("Errors", errors, "\n"); cat("STDbeta", stdbeta, "\n") }
  nom <- 0; denom <- 0; group <- 1
  N <- length(t1)
  J <- length(munique(geno))
  for(g in munique(geno)){
    ind <- which(geno == g)
    Nj <- length(ind)
    t1mean <- mean(t1[ind], na.rm=TRUE)
    value <- 0
    for(x in ind){
      if(!is.na(t1[x])){ value <- value + (t1[x] - t1mean)^2 }
    }
    nom   <- nom + ((betas[group]^2 - stdbeta^2) * value)
    denom <- denom + ((Nj-2)*errors[group])
    group <- group + 1
  }
  if(verbose) cat("SlopeTest Nom/Denom:",nom,"/",denom,"\n")
  nom   <- (1/(J-1)) * nom
  denom <- (1/(N-2*J)) * denom
  res <- NULL
  res$f     <- (nom / denom)
  res$betas <- betas
  return(res)
}

printSignificantSlope <- function(p1, m, p2, pval, slopes, header=FALSE){
  if(!header){
    cat("Phe1\tMar\tPh2\tF\tPvalue", paste0("\tSlope", 1:length(slopes$betas)),"\n")
  }
  cat(p1, "\t", m, "\t", p2, "\t", slopes$f ,"\t", pval, paste0("\t", slopes$betas),"\n")
  return(TRUE);
}

#-- Normal interface --#
scanSlopes <- function(genotypes, phenotypes, pheno.col = 1, doRank = FALSE, verbose = FALSE){
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")

  if(doRank) phenotypes <- apply(phenotypes, 2, rank)
  phenames <- colnames(phenotypes)
  marnames <- colnames(genotypes)
  matrix <- NULL
  header <- FALSE
  for(m in 1:ncol(genotypes)){
    res <- NULL
    for(p in 1:ncol(phenotypes)){
      if(p != pheno.col){ # Do a stdSlopeTest and convert to P
          J       <- length(munique(genotypes[,m]))
          N       <- length(phenotypes[,pheno.col])
          slopes  <- stdSlopeTest(phenotypes[,pheno.col], phenotypes[,p], genotypes[,m], verbose)
          pval    <- 1.0 - pf(slopes$f, J-1, N-2*J)
          if(pval < (0.05 /ncol(genotypes))){ #Print the summary to the screen (TODO: Print to screen / file)
            header <- printSignificantSlope(phenames[pheno.col], marnames[m], phenames[p], pval, slopes, header)
          }
          res <- c(res, pval)
      }else{ res <- c(res, 1.0) } # No chance in hell this will be different for the same trait vs itseld
    }
    matrix <- rbind(matrix, res)
  }
  rownames(matrix) <- marnames
  colnames(matrix) <- phenames
  return(matrix)
}

#-- R/qtl interface --#
scanSlopes.cross <- function(cross, pheno.col = 1, doRank = FALSE, verbose = FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  return(scanSlopes(pull.geno(cross), pull.pheno(cross), pheno.col, doRank, verbose))
}

# end of ctl.slope.R

