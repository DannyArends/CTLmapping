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
    results <- rbind(results,apply(genotypes,2, 
      function(geno){
        linmod <- lm(phenotypes[,x] ~ geno)
        -log10(anova(linmod)[[5]][1])
      }
    ))
    if(verbose){
      cat("Phenotype:",colnames(phenotypes)[x],"\n")
    }
    rnames <- c(rnames,colnames(phenotypes)[x])
    cnt <- cnt +1
  }
  rownames(results) <- rnames
  class(results) <- c(class(results),"QTLscan")
  results
}

# end of qtl.scan.R
