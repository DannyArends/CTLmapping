#
# ctl.scan.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Jan, 2011
# 
# R functions to do CTL mapping
#

#-- CTLscan main function --#
CTLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), method = c("Exact","Power","Adjacency"), n.perms=100, strategy = c("Full", "Pairwise"), conditions = NULL, have.qtls = NULL, n.cores = 1, geno.enc = c(1,2), verbose = FALSE){
  if(missing(genotypes) || is.null(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  st <- proc.time()
  if(verbose) cat("Checking data\n")
  toremove <- check.genotypes(genotypes, geno.enc, verbose)
  phenotypes <- apply(phenotypes,2,rank)
  results <- vector("list",length(pheno.col))
  if(verbose) cat("data.overview:",ncol(phenotypes)," phenotypes,", nrow(phenotypes)," individuals, ", ncol(genotypes)," markers\n")

  if(!is.null(toremove)){
    cat("Cleaning genotype data for mapping\n")
    if(length(toremove) == ncol(genotypes)) stop("Analysis would remove all markers\n")
    warning(paste("Removing",length(toremove),"/",ncol(genotypes),"genotype markers")) 
    genotypes <- genotypes[,-toremove]    
  }
  if(n.cores==1){
    idx <- 1
    for(phe in pheno.col){
      results[[idx]] <- CTLmapping(genotypes, phenotypes, pheno.col=phe, method = method, n.perms = n.perms, strategy=strategy, have.qtls = have.qtls, geno.enc=geno.enc, verbose=verbose)
      idx <- idx + 1
    }
  }else{
    min.cores <- min(n.cores,length(pheno.col))
    if(min.cores != n.cores) warning("Changed the number of cores available (",n.cores," to ",min.cores,") for parallel computation")
    cl <- parallel::makeCluster(rep("localhost", min.cores))
    parallel::clusterEvalQ(cl, library(ctl))
    results <- parallel::parLapply(cl, pheno.col, function(x, have.qtls){
      CTLmapping(genotypes, phenotypes, pheno.col=x, method = method, n.perms = n.perms, strategy=strategy, have.qtls = have.qtls, geno.enc=geno.enc, verbose=verbose)
    },have.qtls)
    parallel::stopCluster(cl)
  }
  if(verbose) cat("Done after: ",(proc.time()-st)[3]," seconds\n")
  class(results) <- c(class(results),"CTLobject")
  results
}

CTLmapping <- function(genotypes, phenotypes, pheno.col=1, method = c("Exact","Power","Adjacency"), n.perms=100, strategy = c("Full", "Pairwise"), have.qtls = NULL, geno.enc=c(1,2), verbose = FALSE){
  if(missing(genotypes) || is.null(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  n.ind = nrow(genotypes); n.mar = ncol(genotypes); n.phe = ncol(phenotypes)
  genotypes[genotypes==geno.enc[1]] <- 0
  genotypes[genotypes==geno.enc[2]] <- 1
  res <- list()
  ss  <- proc.time()
  if(!is.null(have.qtls)){
    res$qtl <- have.qtls[,pheno.col]
  }else{
    tryCatch(res$qtl <- apply(genotypes,2,function(x){-log10(anova(lm(phenotypes[,pheno.col] ~ x))[[5]][1])}), error = function(e) e)
  }
  if(any(is.na(genotypes)))  genotypes[is.na(genotypes)]   <- -999
  if(any(is.na(phenotypes))) phenotypes[is.na(phenotypes)] <- -999
  e1 <- proc.time()

  perm.type = 0
  perms = as.double(rep(0,n.perms))
  if(strategy[1] == "Pairwise"){
    perm.type = 1
    perms = as.double(rep(0,n.perms*n.phe))
  }

  alpha <- 1; gamma <- 1
  if(method[1] != "Exact") gamma <- 2
  if(method[1] == "Adjacency") alpha <- 2

	result <- .C("R_mapctl",as.integer(n.ind), as.integer(n.mar), as.integer(n.phe),
                    			as.integer(unlist(genotypes)), as.double(unlist(phenotypes)),
                          as.integer((pheno.col-1)), as.integer(n.perms),
                          as.integer(alpha),as.integer(gamma),
                          as.integer(perm.type),
                          dcor =as.double(rep(0,n.mar*n.phe)),
                          perms=as.double(perms),
                          ctl  =as.double(rep(0,n.mar*n.phe)),
                          as.integer(verbose), 
                          PACKAGE="ctl")
  e2 <- proc.time()
  res$dcor  <- matrix(result$dcor, n.mar, n.phe)
  res$perms <- result$perms
  if(perm.type==1){
    res$perms <- matrix(result$perms, n.perms, n.phe)
  }
  if(method[1] != "Exact"){ # Likelihoods are computed in C
    res$ctl   <- matrix(result$ctl, n.mar, n.phe)
  }else{  #Exact
    res$ctl   <- (2*dnorm(matrix(result$dcor, n.mar, n.phe))) * (n.mar * n.phe) #Single trait bonferonni correction
    res$ctl[res$ctl > 1] <- 1
    res$ctl   <- -log10(res$ctl)
  }
  if(any(is.na(res$dcor))){
    warning("NaN DCOR scores detected, no variance ?")
    #res$ctl[is.na(res$dcor)] <- 0
  }
  res$ctl  <- Matrix(res$ctl); res$dcor <- Matrix(res$dcor) # Use sparse matrix package to save memory
  rownames(res$dcor) <- colnames(genotypes); colnames(res$dcor) <- colnames(phenotypes)
  rownames(res$ctl)  <- colnames(genotypes); colnames(res$ctl)  <- colnames(phenotypes)
  attr(res,"name") <- colnames(phenotypes)[pheno.col]
  class(res) <- c(class(res),"CTLscan")
  if(verbose) cat("Phenotype ",pheno.col,": Done after ",as.numeric(e2[3]-e1[3])," ",as.numeric(e1[3]-ss[3])," seconds\n",sep="")
  invisible(res)
}

#-- R/qtl interface --#
CTLscan.cross <- function(cross, pheno.col, method = c("Exact","Power","Adjacency"), n.perms=100, strategy = c("Full", "Pairwise"), conditions=NULL, have.qtls=NULL, n.cores=2, verbose=FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(has_rqtl()){
    require(qtl)

    geno.enc <- NULL
    if(any(class(cross)=="bc") || any(class(cross)=="riself") || any(class(cross)=="risib")) geno.enc <- c(1,2)
    if(is.null(geno.enc)) stop("Class of the cross needs to be either: riself,risib or bc")
    rqtl_pheno <- qtl::pull.pheno(cross)
    rqtl_c <- which(colnames(rqtl_pheno) %in% c("sex","pgm")) #R/qtl adds additional columns sex and pgm
    if(length(rqtl_c) > 0) rqtl_pheno <- rqtl_pheno[,-rqtl_c]
    phenotypes <- apply(rqtl_pheno,2,as.numeric)         #R/qtl phenotypes are a data.frame (Need matrix)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    genotypes <- qtl::pull.geno(cross)
    CTLscan(genotypes=genotypes, phenotypes=phenotypes, pheno.col=pheno.col, method=method, n.perms=n.perms, strategy=strategy,
            conditions=conditions, have.qtls=have.qtls, n.cores=n.cores, geno.enc=geno.enc, verbose=verbose)
  }else{
    warning("Please install the R/qtl library (install.packages(\"qtl\"))")
  }
}

# end of ctl.scan.R
