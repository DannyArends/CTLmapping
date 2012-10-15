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
CTLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), n.perms=100, conditions = NULL, have.qtls = NULL, n.cores = 1, geno.enc = c(1,2), verbose = FALSE){
  if(missing(genotypes) || is.null(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  st <- proc.time()
  cat("Stage 0.0: Checking data\n")
  toremove <- check.genotypes(genotypes, geno.enc, verbose)

  results <- vector("list",length(pheno.col))
  stage <- 1
  cat("data.overview:",ncol(phenotypes)," phenotypes,", nrow(phenotypes)," individuals, ", ncol(genotypes)," markers\n")

  if(!is.null(toremove)){
    cat("Stage 0.",stage,": Cleaning genotype data for mapping\n",sep="")
    warning(paste("Removing",length(toremove),"/",ncol(genotypes),"genotype markers")) 
    genotypes <- genotypes[,-toremove]
  }
  idx <- 1
  for(phe in pheno.col){
    qtl <- NULL
    if(!is.null(have.qtls) && idx < ncol(have.qtls)) qtl <- have.qtls[,idx]
    results[[idx]] <- CTLmapping(genotypes, phenotypes, pheno.col=phe, n.perms = n.perms, has.qtl = qtl, geno.enc=geno.enc, verbose=verbose)
    attr(results[[idx]],"name") <- colnames(phenotypes)[phe]
    class(results[[idx]]) <- c(class(results[[idx]]),"CTLscan")
    idx <- idx + 1
  }
  cat("Done after: ",(proc.time()-st)[3]," seconds\n")
  class(results) <- c(class(results),"CTLobject")
  results
}

CTLmapping <- function(genotypes, phenotypes, pheno.col=1, n.perms=100, has.qtl = NULL, geno.enc=c(1,2), verbose = FALSE){
  if(missing(genotypes) || is.null(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  n.ind = nrow(genotypes); n.mar = ncol(genotypes); n.phe = ncol(phenotypes)
  genotypes[genotypes==geno.enc[1]] <- 0
  genotypes[genotypes==geno.enc[2]] <- 1
  res <- list()
  ss  <- proc.time()
  if(!is.null(has.qtl)){
    res$qtl <- has.qtl
  }else{
    res$qtl <- apply(genotypes,2,function(x){-log10(cor.test(x, phenotypes[,pheno.col])$p.value)})
  }
  if(any(is.na(genotypes)))  genotypes[is.na(genotypes)]   <- -999
  if(any(is.na(phenotypes))) phenotypes[is.na(phenotypes)] <- -999
  e1 <- proc.time()
	result <- .C("R_mapctl",as.integer(n.ind), as.integer(n.mar), as.integer(n.phe),
                    			as.integer(unlist(genotypes)), as.double(unlist(phenotypes)),
                          as.integer((pheno.col-1)), as.integer(n.perms),
                          dcor=as.double(rep(0,n.mar*n.phe)),
                          perms=as.double(rep(0,n.perms)),
                          ctl=as.double(rep(0,n.mar*n.phe)), 
                          PACKAGE="ctl")
  e2 <- proc.time()
  res$dcor  <- apply(matrix(result$dcor, n.mar, n.phe), 2, round, digits=2)
  res$perms <- result$perms
  res$ctl   <- Matrix(result$ctl, n.mar, n.phe) # Use matrix package to save memory
  if(any(is.na(res$dcor))){
    warning("NaN DCOR scores detected, no variance ?")
    res$ctl[is.na(res$dcor)] <- 0
  }
  rownames(res$dcor) <- colnames(genotypes); colnames(res$dcor) <- colnames(phenotypes)
  rownames(res$ctl)  <- colnames(genotypes); colnames(res$ctl)  <- colnames(phenotypes)
  attr(res,"name") <- colnames(phenotypes)[pheno.col]
  class(res) <- c(class(res),"CTLscan")
  if(verbose) cat("Phenotype ",pheno.col,": Done after ",as.numeric(e2[3]-e1[3])," ",as.numeric(e1[3]-ss[3])," seconds\n",sep="")
  invisible(res)
}

#-- R/qtl interface --#
CTLscan.cross <- function(cross, pheno.col, n.perms=100, conditions = NULL, have.qtls = NULL, n.cores=2, verbose = FALSE){
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
    CTLscan(genotypes=genotypes, phenotypes=phenotypes, pheno.col=pheno.col, n.perms=n.perms, 
            conditions=conditions, have.qtls=have.qtls, n.cores=n.cores, geno.enc=geno.enc, verbose=verbose)
  }else{
    warning("Please install the R/qtl library (install.packages(\"qtl\"))")
  }
}

# end of ctl.scan.R
