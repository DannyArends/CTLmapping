#
# ctl.scan.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Feb, 2012
# first written Jan, 2011
# 
# R functions to do CTL mapping
#

#-- CTLscan main function --#
CTLscan <- function(genotypes, phenotypes, geno.enc=c(1,2), pheno.col = 1:ncol(phenotypes), have.qtl, method = c("pearson", "kendall", "spearman"), conditions = NULL, n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  
  cat("Stage 0.0: Checking data\n")
  toremove <- check.genotypes(genotypes, geno.enc, verbose)
  results <- vector("list",length(pheno.col))
  stage <- 1
  if(missing(have.qtl)){
    cat("Stage 0.1: Mapping Trait - Marker associations (QTL)\n")
    attr(results,"qtl") <- QTLscan(genotypes, phenotypes, conditions, verbose=verbose)$qtl
    stage <- 2
  }else{
    if(ncol(have.qtl) != ncol(genotypes)) cat("[SEVERE] argument 'have.qtl' should be of size:",ncol(phenotypes),ncol(genotypes))
    if(nrow(have.qtl) != ncol(phenotypes)) cat("[SEVERE] argument 'have.qtl' should be of size:",ncol(phenotypes),ncol(genotypes))
    attr(results,"qtl") <- have.qtl
  }
  if(!is.null(toremove)){
    cat("Stage 0.",stage,": Cleaning genotype data for mapping\n",sep="")
    warning(paste("Removing",length(toremove),"/",ncol(genotypes),"genotype markers")) 
    genotypes <- genotypes[,-toremove]
  }
  idx <- 1
  for(p in pheno.col){
    cat("Stage ",idx,".0: Mapping Correlated Traits Loci (CTL)\n",sep="")
    results[[idx]]$ctl <- CTLmapping(genotypes, phenotypes, geno.enc=geno.enc, p, method=method, verbose)
    results[[idx]]$qtl <- attr(results,"qtl")[p, ]
    if(n.perm > 0){
      cat("Stage ",idx,".1: Permutation\n",sep="")
      results[[idx]]$p <- CTLpermute(genotypes, phenotypes, geno.enc=geno.enc, p, method=method, n.perm, n.cores, directory, saveFiles, verbose)
      
      if(verbose)cat("Stage ",idx,".2: Transformation into LOD\n",sep="")
      results[[idx]]$l <- toLod(results[[idx]], FALSE, verbose)
    }else{
      cat("Stage ",idx,".1: Skipping permutation\n",sep="")
      if(verbose)cat("Stage ",idx,".2: Skipping transformation into LOD\n",sep="")
    }
    class(results[[idx]]) <- c(class(results[[idx]]),"CTLscan")
    idx <- idx + 1
  }
  class(results) <- c(class(results),"CTLobject")
  results
}

CTLmapping <- function(genotypes, phenotypes, geno.enc=c(1,2), pheno.col = 1, method = c("pearson", "kendall", "spearman"), verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  if(length(pheno.col)!=1) stop("CTLmapping can only scan 1 phenotype at once, use CTLscan for multiple phenotypes")
  ss <- proc.time()
  results <- NULL
  ctlprofile <- apply(genotypes,2, 
    function(geno){
      geno1 <- which(geno==geno.enc[1])
      geno2 <- which(geno==geno.enc[2])
      cor1 <- cor(phenotypes[geno1,pheno.col],phenotypes[geno1,],use="pair",method=method[1])
      cor2 <- cor(phenotypes[geno2,pheno.col],phenotypes[geno2,],use="pair",method=method[1])
      #sign(cor1)*(cor1^2)-sign(cor2)*(cor2^2)
      (cor1 - cor2)^2
    }
  )
  rownames(ctlprofile) <- colnames(phenotypes)
  colnames(ctlprofile) <- colnames(genotypes)
  results <- ctlprofile
  attr(results,"name") <- colnames(phenotypes)[pheno.col]
  class(results) <- c(class(results),"CTL")
  ee <- proc.time()
  if(verbose){
    cat("  - CTLscan of",colnames(phenotypes)[pheno.col],"took",as.numeric(ee[3]-ss[3]),"seconds\n")
  }
  results
}

#-- R/qtl interface --#
CTLscan.cross <- function(cross, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(has_rqtl()){
    require(qtl)

    geno.enc <- NULL
    if(any(class(cross)=="bc") || any(class(cross)=="riself") || any(class(cross)=="risib")) geno.enc <- c(1,2)
    if(is.null(geno.enc)) stop("Class of the cross needs to be either: riself,risib or bc")

    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    genotypes <- qtl::pull.geno(cross)
    CTLscan(genotypes = genotypes, phenotypes =  phenotypes, geno.enc=geno.enc, 
            pheno.col = pheno.col, method=method, n.perm=n.perm, n.cores=n.cores, 
            directory=directory, saveFiles=saveFiles, verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}

# end of ctl.scan.R
