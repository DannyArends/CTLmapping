#
# qcl.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Jan, 2011
# 
# R functions to do QCL mapping
#

#-- QCLscan main function --#
QCLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), method = c("pearson", "kendall", "spearman"), n.perm=0, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE){
  results <- vector("list",length(pheno.col))
  idx <- 1
  cat(method,"\n")
  for(p in pheno.col){
    cat("Stage 1: Scanning QCL\n")
    results[[idx]]$qcl <- QCLmapping(genotypes, phenotypes, p, method=method, verbose)
    cat("Stage 2: Scanning QTL\n")
    results[[idx]]$qtl <- QTLscan(genotypes, phenotypes, p, verbose)
    
    if(n.perm > 0){
      cat("Stage 3: Permutation\n")
      results[[idx]]$p <- QCLpermute(genotypes, phenotypes, p, method=method, n.perm, n.cores, directory, saveFiles, verbose)
      
      cat("Stage 4: Transformation into LOD\n")
      results[[idx]]$l <- toLod(results[[idx]], TRUE, FALSE)
    }else{
      cat("Stage 3: Skipping permutation\n")
      cat("Stage 4: Skipping transformation into LOD\n")
    }
    class(results[[idx]]) <- c(class(results[[idx]]),"QCLscan")
    idx <- idx + 1
  }
  class(results) <- c(class(results),"QCLobject")
  results
}

QCLmapping <- function(genotypes, phenotypes, pheno.col = 1, method = c("pearson", "kendall", "spearman"), verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  ss <- proc.time()
  results <- NULL
  profile <- apply(genotypes,2, 
    function(geno){
      cor1 <- cor(phenotypes[geno==1,pheno.col],phenotypes[geno==1,],use="pair",method=method[1])
      cor2 <- cor(phenotypes[geno==2,pheno.col],phenotypes[geno==2,],use="pair",method=method[1])
      sign(cor1)*(cor1^2)-sign(cor2)*(cor2^2)
    }
  )
  rownames(profile) <- colnames(phenotypes)
  colnames(profile) <- colnames(genotypes)
  results <- profile
  attr(results,"name") <- colnames(phenotypes)[pheno.col]
  class(results) <- c(class(results),"QCL")
  ee <- proc.time()
  if(verbose){
    cat("  - QCLscan of",colnames(phenotypes)[pheno.col],"took",as.numeric(ee[3]-ss[3]),"seconds\n")
  }
  results
}

#-- R/qtl interface --#
QCLscan.cross <- function(cross, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(has_rqtl()){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    genotypes <- qtl::pull.geno(cross)
    QCLscan(genotypes, phenotypes, pheno.col, method, n.perm, n.cores, directory, saveFiles, verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}
