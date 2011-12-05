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
QCLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE){
  results <- vector("list",length(pheno.col))
  idx <- 1
  for(p in pheno.col){
    cat("Scanning\n")
    results[[idx]]$s <- QCLmapping(genotypes, phenotypes, p, verbose)
    cat("Permutation\n")
    results[[idx]]$p <- QCLpermute(genotypes, phenotypes, p, n.perm, n.cores, directory, saveFiles, verbose)
    cat("toLOD\n")
    results[[idx]]$l <- toLod(results[[idx]], FALSE)
    class(results[[idx]]) <- c(class(results[[idx]]),"QCLscan")
    idx <- idx + 1
  }
  class(results) <- c(class(results),"QCLobject")
  results
}

QCLmapping <- function(genotypes, phenotypes, pheno.col = 1, verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  results <- NULL
  profile <- apply(genotypes,2, 
    function(geno){
      cor1 <- cor(phenotypes[geno==1,pheno.col],phenotypes[geno==1,],use="pair")
      cor2 <- cor(phenotypes[geno==2,pheno.col],phenotypes[geno==2,],use="pair")
      sign(cor1)*(cor1^2)-sign(cor2)*(cor2^2)
    }
  )
  rownames(profile) <- colnames(phenotypes)
  colnames(profile) <- colnames(genotypes)
  results <- profile
  attr(results,"name") <- colnames(phenotypes)[pheno.col]
  class(results) <- c(class(results),"QCL")
  if(verbose){
    cat("Phenotype:",colnames(phenotypes)[pheno.col],"\n")
  }
  results
}

#-- R/qtl interface --#
QCLscan.cross <- function(cross, pheno.col, verbose = FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(get(".has_rqtl", envir = .QclEnv)){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    genotypes <- qtl::pull.geno(cross)
    QCLscan(genotypes, phenotypes, pheno.col, verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}
