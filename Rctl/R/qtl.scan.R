#
# qtl.scan.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Jan, 2011
# 
# R functions to do QTL mapping
#

#-- QTLscan main function --#

QTLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), verbose = FALSE){
  if(missing(genotypes)) stop("genotypes are missing")
  if(missing(phenotypes)) stop("phenotypes are missing")
  results <- NULL
  rnames <- NULL
  cnt <- 1
  for(x in pheno.col){
    ss <- proc.time()
    results <- rbind(results,apply(genotypes,2, 
      function(geno){
        linmod <- lm(phenotypes[,x] ~ geno)
        -log10(anova(linmod)[[5]][1])
      }
    ))
    ee <- proc.time()
    if(verbose){
      cat("  - QTLscan of",colnames(phenotypes)[x],"took",as.numeric(ee[3]-ss[3]),"seconds\n")
    }
    rnames <- c(rnames,colnames(phenotypes)[x])
    cnt <- cnt +1
  }
  rownames(results) <- rnames
  class(results) <- c(class(results),"QTLscan")
  results
}


#-- R/qtl interface --#
QTLscan.cross <- function(cross, pheno.col, verbose = FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(has_rqtl()){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    genotypes <- qtl::pull.geno(cross)
    QTLscan(genotypes, phenotypes, pheno.col, verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}

# end of qtl.scan.R
