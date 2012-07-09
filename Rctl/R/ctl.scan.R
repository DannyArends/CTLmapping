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
CTLscan <- function(genotypes, phenotypes, geno.enc=c(1,2), pheno.col = 1:ncol(phenotypes), have.qtl, method = c("pearson","pearsonordered", "kendall", "spearman"), conditions = NULL, n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  
  cat("Stage 0.0: Checking data\n")
  toremove <- check.genotypes(genotypes, geno.enc, verbose)
  results <- vector("list",length(pheno.col))
  stage <- 1
  cat("[] ",nrow(phenotypes)," individuals, ",ncol(genotypes)," markers\n")
  if(missing(have.qtl)){
    cat("Stage 0.1: Mapping Trait - Marker associations (QTL)\n")
    attr(results,"qtl") <- QTLscan(genotypes, phenotypes, pheno.col, conditions, n.cores, verbose=verbose)$qtl
    stage <- 2
  }else{
    if(ncol(have.qtl) != ncol(genotypes))   cat("[SEVERE] argument 'have.qtl' should be of size:",length(pheno.col)," ",ncol(genotypes),"\n")
    if(nrow(have.qtl) != length(pheno.col)) cat("[SEVERE] argument 'have.qtl' should be of size:",length(pheno.col)," ",ncol(genotypes),"\n")
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
      res <- matrix(0,1,ncol(phenotypes))
      for(x in 1:(length(geno.enc)-1)){
        for(y in (x+1):length(geno.enc)){
          if(x != y){
            geno1 <- which(geno==geno.enc[x])
            geno2 <- which(geno==geno.enc[y])
            cor1 <- cor(phenotypes[geno1,pheno.col],phenotypes[geno1,],use="pair",method=method[1])
            cor2 <- cor(phenotypes[geno2,pheno.col],phenotypes[geno2,],use="pair",method=method[1])
            res = res + (cor1 - cor2)^2
          }
        }
      }
      return(res)
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
CTLscan.cross <- function(cross, pheno.col, have.qtl, method = c("pearson", "pearsonordered", "kendall", "spearman"), n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE){
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
    CTLscan(genotypes = genotypes, phenotypes =  phenotypes, geno.enc=geno.enc, 
            pheno.col = pheno.col, have.qtl=have.qtl, method=method, n.perm=n.perm, n.cores=n.cores, 
            directory=directory, saveFiles=saveFiles, verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}

# end of ctl.scan.R
