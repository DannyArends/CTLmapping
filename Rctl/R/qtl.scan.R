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
  qtldata <- lapply(models,"[",what)
  res <- matrix(unlist(qtldata),length(phenonames),length(genonames))
  rownames(res) <- phenonames
  colnames(res) <- genonames
  if(do.log) res <- -log10(res)
  res
}

#-- clean.phenotypes function to prepare phenotypes for QTL mapping--#
clean.phenotypes <- function(phenotypes, verbose = TRUE){
  phenotypes <- as.matrix(apply(phenotypes,2,as.numeric))
  torem <- which(apply(apply(phenotypes,2,is.na),2,sum)==nrow(phenotypes))
  if(!is.na(torem&&1)){
    if(verbose) warning("Filling ",length(torem)," phenotypes (",paste(torem,collapse=","),") that are NA with 1s\n")
    for(x in torem){ phenotypes[,x] <- rep(1,nrow(phenotypes)) }
  }
  if(is.null(colnames(phenotypes))) colnames(phenotypes) <- paste("Pheno",1:ncol(phenotypes),sep="")
  phenotypes
}

#-- QTLscan main function --#
map.fast <- function(marker, phenotypes, conditions = NULL, verbose = FALSE){
  st  <- proc.time()
  res <- NULL
  if(is.null(conditions)){
    models  <- aov(phenotypes ~ marker)
    modelinfo <- summary(models)
    res$qtl <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
    res$eff <- unlist(models$coefficients[2,])
  }else{
    models  <- aov(phenotypes ~ conditions + marker + conditions:marker)
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
  phenotypes <- clean.phenotypes(phenotypes)
  if("snow" %in% rownames(installed.packages()) && n.core > 1){
    require("snow")
    cl <- snow::makeCluster(n.core)
    models <- snow::parApply(cl, genotypes, 2, "map.fast", phenotypes, conditions)
    snow::stopCluster(cl)
  }else{
    models <- apply(genotypes, 2, "map.fast", phenotypes, conditions)  
  }
  cat("Stage 0.1: Done with QTL mapping\n")
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
