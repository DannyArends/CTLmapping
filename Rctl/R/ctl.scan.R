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
CTLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), method = c("pearson", "kendall", "spearman"), n.perm=100, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  
  cat("Stage 0.0: Checking data\n")
  toremove <- check.genotypes(genotypes, genotype.values, verbose)
  cat("Stage 0.1: Mapping Trait - Marker associations (QTL)\n")
  results <- vector("list",length(pheno.col))
  attr(results,"qtl") <- QTLscan(genotypes, phenotypes, verbose=verbose)$qtl
  cat("Stage 0.2: Cleaning data for CTL mapping\n")
  if(!is.null(toremove)) genotypes <- genotypes[,-toremove]

  idx <- 1
  for(p in pheno.col){
    cat("Stage ",idx,".0: Mapping Correlated Traits Loci - CTL\n",sep="")
    results[[idx]]$ctl <- CTLmapping(genotypes, phenotypes, p, method=method, genotype.values, verbose)
    results[[idx]]$qtl <- attr(results,"qtl")[p, ]
    if(n.perm > 0){
      cat("Stage ",idx,".1: Permutation\n",sep="")
      results[[idx]]$p <- CTLpermute(genotypes, phenotypes, p, method=method, n.perm, n.cores, genotype.values, directory, saveFiles, verbose)
      
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

CTLmapping <- function(genotypes, phenotypes, pheno.col = 1, method = c("pearson", "kendall", "spearman"), genotype.values=c(1,2), verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  if(length(pheno.col)!=1) stop("CTLmapping can only scan 1 phenotype at once, use CTLscan for multiple phenotypes")
  ss <- proc.time()
  results <- NULL
  ctlprofile <- apply(genotypes,2, 
    function(geno){
      geno1 <- which(geno==genotype.values[1])
      geno2 <- which(geno==genotype.values[2])
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
CTLscan.cross <- function(cross, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=100, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose = FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(has_rqtl()){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    genotypes <- qtl::pull.geno(cross)
    CTLscan(genotypes = genotypes, phenotypes =  phenotypes, pheno.col = pheno.col, 
            method=method, n.perm=n.perm, n.cores=n.cores,genotype.values=genotype.values, 
            directory=directory, saveFiles=saveFiles, verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}

# end of ctl.scan.R
