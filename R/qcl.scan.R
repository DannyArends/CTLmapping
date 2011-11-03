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
QCLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), verbose = FALSE){
  if(missing(genotypes)) stop("genotypes are missing")
  if(missing(phenotypes)) stop("phenotypes are missing")
  results <- vector("list",length(pheno.col))
  cnt <- 1
  for(x in pheno.col){
    profile <- apply(genotypes,2, 
      function(geno){
        cor1 <- cor(phenotypes[geno==1,x],phenotypes[geno==1,],use="pair")
        cor2 <- cor(phenotypes[geno==2,x],phenotypes[geno==2,],use="pair")
        sign(cor1)*(cor1^2)-sign(cor2)*(cor2^2)
      }
    )
    rownames(profile) <- colnames(phenotypes)
    colnames(profile) <- colnames(genotypes)
    results[[cnt]] <- profile
    attr(results[[cnt]],"name") <- colnames(phenotypes)[x]
    class(results[[cnt]]) <- c(class(results[[cnt]]),"QCL")
    if(verbose){
      cat("Phenotype:",colnames(phenotypes)[x],"\n")
    }
    cnt <- cnt +1
  }
  class(results) <- c(class(results),"QCLscan")
  results
}

#-- R/qtl interface --#
QCLscan.cross <- function(cross, pheno.col, verbose = FALSE){
  if(missing(cross)) stop("cross is missing")
  if(get(".has_rqtl", envir = .QclEnv)){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    genotypes <- qtl::pull.geno(cross)
    QCLscan(genotypes, phenotypes,pheno.col=pheno.col,verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}
