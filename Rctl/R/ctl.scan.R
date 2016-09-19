#
# ctl.scan.R
#
# copyright (c) 2010-2013 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Feb, 2013
# first written Jan, 2011
# 
# R functions to do CTL mapping
#

#-- CTLscan main function --#
CTLscan <- function(genotypes, phenotypes, phenocol, nperm = 100, strategy = c("Exact", "Full", "Pairwise"), conditions = NULL, qtls = NULL, ncores = 1, parametric = FALSE, verbose = FALSE){
  st <- proc.time()
  if(missing(genotypes) || is.null(genotypes))  stop("argument 'genotypes' is missing, with no default")
  if(any(class(genotypes)=="cross")){
    CTLscan.cross(genotypes, phenocol=phenocol, nperm=nperm, strategy=strategy, conditions=conditions,
                  qtls=qtls, ncores=ncores, parametric = parametric, verbose=verbose)
  }
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  if(missing(phenocol)) phenocol = 1:ncol(phenotypes)
  n.phe      <- ncol(phenotypes); n.ind1 <- nrow(phenotypes)
  n.mar      <- ncol(genotypes);  n.ind2 <- nrow(genotypes)
  if(!is.null(qtls) && ncol(qtls) != length(phenocol)) stop("Number of QTLs doesn't match")
  if(n.ind1 != n.ind2) stop("Number of individuals doesn't match between genotypes & phenotypes")
  if(n.phe < 2) stop("argument 'phenotypes' needs at least 2 columns")
  #if(length(phenocol) > n.phe) stop(paste0("trying to scan more phenotypes then available (", length(phenocol), "/", n.phe, ")"))
  if(verbose){
    cat("Data", n.phe ,"phenotypes,", paste0(n.ind1, "/", n.ind2), "individuals,", n.mar ,"markers\n")
    cat("Data checks finished after:",(proc.time()-st)[3],"seconds\n")
  }
  if(!parametric){
    phenotypes <- apply(phenotypes, 2, rank) # Parametric vs Non-Parametric testing
    if(verbose) cat("Data ranking finished after:",(proc.time()-st)[3],"seconds\n")
  }
  results    <- vector("list",length(phenocol))

  if(ncores == 1){ # Single core debug call
    idx <- 1
    for(phe in phenocol){
      results[[idx]] <- CTLmapping(genotypes, phenotypes, phenocol=phe, nperm=nperm,
                             strategy=strategy, qtls=qtls, verbose=verbose)
      idx <- idx + 1
    }
  }else{ # Normally we would want to use multiple cores
    min.cores <- min(ncores, length(phenocol))
    if(min.cores != ncores) warning("Reduced n.cores (", ncores, "to", min.cores, ")")
    cl <- parallel::makeCluster(rep("localhost", min.cores))
    parallel::clusterEvalQ(cl, library(ctl))
    results <- parallel::parLapply(cl, phenocol, function(x, qtls){
      CTLmapping(genotypes, phenotypes, phenocol=x, nperm=nperm, strategy=strategy,
                 qtls = qtls, verbose=verbose)
    }, qtls)
    parallel::stopCluster(cl)
  }
  if(verbose) cat("Done after:",(proc.time()-st)[3],"seconds\n")
  class(results) <- c(class(results),"CTLobject")
  invisible(results)
}

CTLmapping <- function(genotypes, phenotypes, phenocol = 1, nperm = 100, strategy = c("Exact", "Full", "Pairwise"), qtls = NULL, verbose = FALSE){
  if(missing(genotypes) || is.null(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")

  n.ind = nrow(genotypes); n.mar = ncol(genotypes); n.phe = ncol(phenotypes)
  res <- list()
  ss  <- proc.time()
  if(!is.null(qtls)){
    res$qtl <- qtls[,phenocol]
  }else{
    res$qtl <- rep(0,ncol(genotypes))
    tryCatch( #Maps QTL profile using a relatively 'slow' approach
      res$qtl <- apply(genotypes, 2, function(x){
                   -log10(anova(lm(phenotypes[,phenocol] ~ x))[[5]][1])
                 }), error = function(e) e)
  }

  #In C we use -999 for missing genotype data
  if(any(is.na(genotypes)))  genotypes[is.na(genotypes)]   <- -999
  if(any(is.na(phenotypes))) phenotypes[is.na(phenotypes)] <- -999
  e1 <- proc.time()

  #Setup return structures for the permutations
  perm.type <- 0
  perms     <- as.double(rep(0, nperm))
  if(strategy[1] == "Full"){ perm.type <- 1; }
  if(strategy[1] == "Pairwise"){
    perm.type = 2; perms = as.double(rep(0,nperm * n.phe))
  }

  result <- .C("R_mapctl",as.integer(n.ind), as.integer(n.mar), as.integer(n.phe),
                          as.integer(unlist(genotypes)), as.double(unlist(phenotypes)),
                          as.integer((phenocol-1)), as.integer(nperm),
                          as.integer(perm.type),
                          dcor =as.double(rep(0,n.mar*n.phe)),
                          perms=as.double(perms),
                          ctl  =as.double(rep(0,n.mar*n.phe)),
                          as.integer(verbose), 
                          PACKAGE="ctl")
  e2 <- proc.time()
  res$dcor  <- matrix(result$dcor, n.mar, n.phe)   # Store the DCOR/Z scores
  res$perms <- result$perms                        # Permutations
  if(perm.type == 1){ res$perms <- matrix(result$perms, nperm, n.phe); }

  res$ctl   <- matrix(result$ctl, n.mar, n.phe)    # CTL likelihoods
  if(any(is.na(res$dcor))) warning("Differential correlation: NaN scores, no variance ?")
  rownames(res$dcor) <- colnames(genotypes); colnames(res$dcor) <- colnames(phenotypes)
  rownames(res$ctl)  <- colnames(genotypes); colnames(res$ctl)  <- colnames(phenotypes)
  attr(res,"name") <- colnames(phenotypes)[phenocol]
  class(res) <- c(class(res),"CTLscan")
  if(verbose){
    t1 <- as.numeric(e2[3]-e1[3]) ; t2 <- as.numeric(e1[3]-ss[3])
    cat("Phenotype ",phenocol,": Done after ", t1," ", t2," seconds\n", sep="")
  }
  invisible(res)
}

#-- R/qtl interface --#
CTLscan.cross <- function(cross, ...){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  rqtl_pheno <- qtl::pull.pheno(cross)
  rqtl_c     <- NULL;               # R/qtl adds additional columns sex and pgm
  cond_id    <- which(colnames(rqtl_pheno) %in% c("sex", "pgm"))
  if(length(cond_id) > 0){
    rqtl_c     <- rqtl_pheno[, cond_id]                      # Add them as conditions
    rqtl_pheno <- rqtl_pheno[,-cond_id]                      # Remove them as phenotypes
  }
  if(ncol(rqtl_pheno) <= 1) stop("Not enough phenotypes in cross object (",ncol(rqtl_pheno), "< 2)")
  phenotypes <- apply(rqtl_pheno, 2, as.numeric)             # R/qtl phenotypes data.frame (need matrix)
  genotypes  <- pull.geno(cross)
  CTLscan(genotypes=genotypes, phenotypes=phenotypes, conditions = rqtl_c, ...)
}

# end of ctl.scan.R

