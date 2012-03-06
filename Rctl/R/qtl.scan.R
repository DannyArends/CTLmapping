#
# qtl.scan.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Jan, 2011
# 
# R functions to do QTL mapping
#

#Recreate a matrix from the per marker modeling results
create.matrix <- function(models, what="qtl", phenonames, genonames, do.log=TRUE){
  res <- matrix(unlist((lapply(models,"[",what))),length(phenonames),length(genonames))
  rownames(res) <- phenonames
  colnames(res) <- genonames
  if(do.log) res <- -log10(res)
  res
}

#-- QTLscan main function --#
map.fast <- function(marker, phenotypes, conditions, verbose = FALSE){
  st  <- proc.time()
  res <- NULL
  if(is.null(conditions)){
    models  <- aov(as.matrix(phenotypes) ~ marker)
    modelinfo <- summary(models)
    res$qtl <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
    res$eff <- unlist(models$coefficients[2,])
  }else{
    models  <- aov(as.matrix(phenotypes) ~ conditions + marker + conditions:marker)
    modelinfo <- summary(models)
    res$env <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
    res$qtl <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
    res$int <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
    res$eff <- unlist(models$coefficients[2,])
  }
  if(verbose) cat("Done in:",(proc.time()-st)[3],"seconds\n")
  res
}

#-- QTLscan main function --#
QTLscan <- function(genotypes, phenotypes, conditions=NULL, n.core=2, verbose = TRUE){
  if(missing(genotypes)) stop("genotypes are missing")
  if(missing(phenotypes)) stop("phenotypes are missing")
  st  <- proc.time()
  if("snow" %in% rownames(installed.packages())){
    require("snow")
    cl <- snow::makeCluster(n.core)
    models <- snow::parApply(cl, genotypes, 2, "map.fast", phenotypes, conditions)
    snow::stopCluster(cl)
  }else{
    models <- apply(genotypes, 2, "map.fast", phenotypes, conditions)  
  }
  res <- NULL
  res$qtl <- create.matrix(models, "qtl", colnames(phenotypes), colnames(genotypes))
  res$eff <- create.matrix(models, "eff", colnames(phenotypes), colnames(genotypes), FALSE)
  if(!is.null(conditions)){
    res$env <- create.matrix(models, "env", colnames(phenotypes), colnames(genotypes))
    res$int <- create.matrix(models, "int", colnames(phenotypes), colnames(genotypes))
  }
  if(verbose)cat("Done in:",(proc.time()-st)[3],"seconds\n")
  res
  class(res) <- c(class(res),"QTLscan")
  res
}

#-- R/qtl interface --#
QTLscan.cross <- function(cross, pheno.col, conditions=NULL, n.core=2, verbose = FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(has_rqtl()){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    genotypes <- qtl::pull.geno(cross)
    QTLscan(genotypes, phenotypes, conditions, n.core, verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}

# end of qtl.scan.R
