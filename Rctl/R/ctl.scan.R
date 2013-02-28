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
CTLscan <- function(genotypes, phenotypes, pheno.col, n.perms = 100, strategy = c("Exact", "Full", "Pairwise"), conditions = NULL, have.qtls = NULL, n.cores = 1, geno.enc = c(1,2), verbose = FALSE){
  if(missing(genotypes) || is.null(genotypes))  stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  if(missing(pheno.col)) pheno.col = 1:ncol(phenotypes)
  st <- proc.time()
  toremove   <- check.genotypes(genotypes, geno.enc, verbose=verbose)
  phenotypes <- apply(phenotypes, 2, rank) # Always use non-parametric statistics
  results    <- vector("list",length(pheno.col))
  n.phe      <- ncol(phenotypes); n.ind1 <- nrow(phenotypes)
  n.mar      <- ncol(genotypes);  n.ind2 <- nrow(genotypes)
  if(verbose){
    cat("Data: phenotypes:", n.phe ," phenotypes, ", n.ind1, " individuals\n")
    cat("Data: genotypes:", n.mar ," markers, ", n.ind2, " individuals\n")
  }
  if(!is.null(have.qtls) && ncol(have.qtls) != length(pheno.col)) stop("Number of QTLs doesn't match")
  if(n.ind1 != n.ind2) stop("Number of individuals doesn't match")
  if(!is.null(toremove)){
    if(length(toremove) == n.mar) stop("Analysis would remove all markers\n")
    warning(paste("Removing genotype markers (", length(toremove),"/", n.mar, ")")) 
    genotypes <- genotypes[,-toremove]
    n.mar     <- ncol(genotypes); 
  }
  if(n.cores==1){
    idx <- 1
    for(phe in pheno.col){
      results[[idx]] <- CTLmapping(genotypes, phenotypes, pheno.col=phe, n.perms=n.perms,
                             strategy=strategy, have.qtls=have.qtls, geno.enc=geno.enc, verbose=verbose)
      idx <- idx + 1
    }
  }else{
    min.cores <- min(n.cores,length(pheno.col))
    if(min.cores != n.cores) warning("Reduced n.cores (",n.cores," to ",min.cores,")")
    cl <- parallel::makeCluster(rep("localhost", min.cores))
    parallel::clusterEvalQ(cl, library(ctl))
    results <- parallel::parLapply(cl, pheno.col, function(x, have.qtls){
      CTLmapping(genotypes, phenotypes, pheno.col=x, n.perms=n.perms, strategy=strategy,
                 have.qtls = have.qtls, geno.enc=geno.enc, verbose=verbose)
    },have.qtls)
    parallel::stopCluster(cl)
  }
  if(verbose) cat("Done after: ",(proc.time()-st)[3]," seconds\n")
  class(results) <- c(class(results),"CTLobject")
  invisible(results)
}

CTLmapping <- function(genotypes, phenotypes, pheno.col = 1, n.perms = 100, strategy = c("Exact", "Full", "Pairwise"), have.qtls = NULL, geno.enc = c(1,2), verbose = FALSE){
  if(missing(genotypes) || is.null(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")

  n.ind = nrow(genotypes); n.mar = ncol(genotypes); n.phe = ncol(phenotypes)
  res <- list()
  ss  <- proc.time()
  if(!is.null(have.qtls)){
    res$qtl <- have.qtls[,pheno.col]
  }else{
    res$qtl <- rep(0,ncol(genotypes))
    #tryCatch( #Maps QTL profile using a 'slow' approach
    #  res$qtl <- apply(genotypes, 2, function(x){
    #               -log10(anova(lm(phenotypes[,pheno.col] ~ x))[[5]][1])
    #             }), error = function(e) e)
  }

  #In C we use -999 for missing genotype data
  if(any(is.na(genotypes)))  genotypes[is.na(genotypes)]   <- -999
  if(any(is.na(phenotypes))) phenotypes[is.na(phenotypes)] <- -999
  e1 <- proc.time()

  #Setup return structures for the permutations
  perm.type = 0
  perms = as.double(rep(0, n.perms))
  if(strategy[1] == "Full"){
    perm.type = 1
  }
  if(strategy[1] == "Pairwise"){
    perm.type = 2
    perms = as.double(rep(0,n.perms * n.phe))
  }

	result <- .C("R_mapctl",as.integer(n.ind), as.integer(n.mar), as.integer(n.phe),
                    			as.integer(unlist(genotypes)), as.double(unlist(phenotypes)),
                          as.integer((pheno.col-1)), as.integer(n.perms),
                          as.integer(1),as.integer(1),  #Setup ALPHA & GAMMA
                          as.integer(perm.type),
                          dcor =as.double(rep(0,n.mar*n.phe)),
                          perms=as.double(perms),
                          ctl  =as.double(rep(0,n.mar*n.phe)),
                          as.integer(verbose), 
                          PACKAGE="ctl")
  e2 <- proc.time()
  #Store the DCOR/Z scores, Permutations, CTLs
  res$dcor  <- matrix(result$dcor, n.mar, n.phe)
  res$perms <- result$perms
  if(perm.type == 1){
    res$perms <- matrix(result$perms, n.perms, n.phe)
  }
  res$ctl   <- matrix(result$ctl, n.mar, n.phe)
  if(any(is.na(res$dcor))) warning("Differential correlation: NaN scores, no variance ?")
  rownames(res$dcor) <- colnames(genotypes); colnames(res$dcor) <- colnames(phenotypes)
  rownames(res$ctl)  <- colnames(genotypes); colnames(res$ctl)  <- colnames(phenotypes)
  attr(res,"name") <- colnames(phenotypes)[pheno.col]
  class(res) <- c(class(res),"CTLscan")
  if(verbose){
     t1 <- as.numeric(e2[3]-e1[3])
     t2 <- as.numeric(e1[3]-ss[3])
     cat("Phenotype ",pheno.col,": Done after ", t1," ", t2," seconds\n", sep="")
  }
  invisible(res)
}

#-- R/qtl interface --#
CTLscan.cross <- function(cross, ...){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  geno.enc <- NULL
  if(any(class(cross)=="bc") || any(class(cross)=="riself") || any(class(cross)=="risib")){
    geno.enc <- c(1,2)
  }else{
    stop("class of cross needs to be either: riself, risib or bc")
  }
  rqtl_pheno <- qtl::pull.pheno(cross)
  rqtl_c     <- NULL;               # R/qtl adds additional columns sex and pgm
  cond_id    <- which(colnames(rqtl_pheno) %in% c("sex", "pgm"))
  if(length(cond_id) > 0){
    rqtl_c     <- rqtl_pheno[, cond_id]                      # Add them as conditions
    rqtl_pheno <- rqtl_pheno[,-cond_id]                      # Remove them as phenotypes
  }  
  phenotypes <- apply(rqtl_pheno, 2, as.numeric)             # R/qtl phenotypes data.frame (need matrix)
  genotypes  <- pull.geno(cross)
  CTLscan(genotypes=genotypes, phenotypes=phenotypes, geno.enc = geno.enc, conditions = rqtl_c, ...)
}

# end of ctl.scan.R
