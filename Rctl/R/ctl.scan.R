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
CTLscan <- function(genotypes, phenotypes, phenocol, nperm = 100, 
                    strategy = c("Exact", "Full", "Pairwise"), ncores = 1, 
                    parametric = FALSE, qtl = TRUE, verbose = FALSE) {

  if(missing(genotypes) || is.null(genotypes))  stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  if(missing(phenocol)) phenocol <- 1:ncol(phenotypes)

  st <- proc.time()

  n.phe <- ncol(phenotypes); n.ind1 <- nrow(phenotypes)
  n.mar <- ncol(genotypes); n.ind2 <- nrow(genotypes)

  if (n.ind1 != n.ind2) {
    stop("number of individuals doesn't match between genotypes (", n.ind2, ") & phenotypes (", n.ind1, ")")
  }
  if (n.phe < 2) {
    stop("argument 'phenotypes' needs at least 2 columns")
  }
  if (length(phenocol) > n.phe) {
    warning(paste0("trying to scan more phenotypes then available (", length(phenocol), ">", n.phe, ")"))
    phenocol <- 1:ncol(phenotypes)
  }
  if (verbose) {
    cat("Data", n.phe ,"phenotypes,", paste0(n.ind1, "/", n.ind2), "individuals,", n.mar ,"markers\n")
    cat("Data checks finished after:",(proc.time()-st)[3],"seconds\n")
  }
  if (!parametric) {
    phenotypes <- apply(phenotypes, 2, rank) # Parametric vs Non-Parametric testing
    if(verbose) cat("Data ranking finished after:",(proc.time()-st)[3],"seconds\n")
  }
  ctlobject <- vector("list", length(phenocol))

  if(ncores == 1){ # Single core debug call
    idx <- 1
    for(phe in phenocol) {
      ctlobject[[idx]] <- CTLmapping(genotypes, phenotypes, phenocol = phe, nperm = nperm,
                                     strategy = strategy, qtl = qtl, verbose = verbose)
      idx <- idx + 1
    }
  }else{ 
    # This is more or less the worst place to do multicore, we duplicate all memory allocated to all cores
    # Normally we would want to use multiple cores, but in the C code, splitting by e.g. marker
    min.cores <- min(ncores, length(phenocol))
    if(min.cores != ncores) warning("Reduced n.cores (", ncores, "to", min.cores, ")")
    cl <- parallel::makeCluster(rep("localhost", min.cores))
    parallel::clusterEvalQ(cl, library(ctl))
    ctlobject <- parallel::parLapply(cl, phenocol, function(x) {
      CTLmapping(genotypes, phenotypes, phenocol = x, nperm = nperm, 
                 strategy = strategy, qtl = qtl, verbose = verbose)
    })
    parallel::stopCluster(cl)
  }
  if(verbose) cat("Done after:", (proc.time()-st)[3], "seconds\n")

  class(ctlobject) <- c(class(ctlobject), "CTLobject")
  invisible(ctlobject)
}

CTLmapping <- function(genotypes, phenotypes, phenocol = 1, nperm = 100, strategy = c("Exact", "Full", "Pairwise"), qtl = TRUE, verbose = FALSE){
  if(missing(genotypes) || is.null(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  ss  <- proc.time()

  n.ind = nrow(genotypes); n.mar = ncol(genotypes); n.phe = ncol(phenotypes)
  ctlscan <- list()
  ctlscan$qtl <- rep(0, n.mar)
  if(qtl) ctlscan$qtl <- QTLmapping(genotypes = genotypes, phenotypes = phenotypes, phenocol = phenocol)

  # In C we use -999 for missing genotype data
  if(any(is.na(genotypes)))  genotypes[is.na(genotypes)]   <- -999
  if(any(is.na(phenotypes))) phenotypes[is.na(phenotypes)] <- -999
  e1 <- proc.time()

  # Setup return structures for the permutations
  perm.type <- 0
  perms <- as.double(rep(0, nperm))
  if(strategy[1] == "Full"){ perm.type <- 1; }
  if(strategy[1] == "Pairwise"){
    perm.type = 2; perms = as.double(rep(0, nperm * n.phe))
  }

  result <- .C("R_mapctl",as.integer(n.ind), as.integer(n.mar), as.integer(n.phe),
                          as.integer(unlist(genotypes)), as.double(unlist(phenotypes)),
                          as.integer((phenocol-1)), as.integer(nperm),
                          as.integer(perm.type),
                          dcor = as.double(rep(0, n.mar * n.phe)),
                          perms = as.double(perms),
                          ctl  = as.double(rep(0, n.mar * n.phe)),
                          as.integer(verbose), 
                          PACKAGE = "ctl")
  e2 <- proc.time()
  ctlscan$dcor  <- matrix(result$dcor, n.mar, n.phe)   # Store the DCOR/Z scores
  ctlscan$perms <- result$perms                        # Permutations
  if(perm.type == 1) { 
    ctlscan$perms <- matrix(result$perms, nperm, n.phe);
  }
  ctlscan$ctl   <- matrix(result$ctl, n.mar, n.phe)    # CTL likelihoods
  if(any(is.na(ctlscan$dcor))) warning("Differential correlation: NaN scores, no variance ?")
  rownames(ctlscan$dcor) <- colnames(genotypes)
  colnames(ctlscan$dcor) <- colnames(phenotypes)

  rownames(ctlscan$ctl)  <- colnames(genotypes)
  colnames(ctlscan$ctl)  <- colnames(phenotypes)

  ctlscan$phenotype <- colnames(phenotypes)[phenocol]     # Set the phenotype name as separate list item
  attr(ctlscan,"name") <- colnames(phenotypes)[phenocol]  # Set it also as attribute to the list

  class(ctlscan) <- c(class(ctlscan),"CTLscan")           # Class is a CTLscan object
  if(verbose){
    t1 <- as.numeric(e2[3]-e1[3]) ; t2 <- as.numeric(e1[3]-ss[3])
    cat("Phenotype ", ctlscan$phenotype, ": Done after ", t1," ", t2, " seconds\n", sep="")
  }
  invisible(ctlscan)
}

#-- R/qtl interface --#
CTLscan.cross <- function(cross, ...){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  rqtl_pheno <- qtl::pull.pheno(cross)
  rqtl_c     <- NULL;               # R/qtl adds additional columns sex and pgm
  cond_id    <- which(colnames(rqtl_pheno) %in% c("sex", "pgm"))
  if(length(cond_id) > 0){
    rqtl_pheno <- rqtl_pheno[,-cond_id]                      # Remove them as phenotypes
  }
  if(ncol(rqtl_pheno) <= 1) stop("Not enough phenotypes in cross object (",ncol(rqtl_pheno), "< 2)")
  phenotypes <- apply(rqtl_pheno, 2, as.numeric)             # R/qtl phenotypes data.frame (need matrix)
  genotypes  <- pull.geno(cross)
  CTLscan(genotypes=genotypes, phenotypes=phenotypes, ...)
}

# end of ctl.scan.R

