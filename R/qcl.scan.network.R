#
# qcl.scan.network.R
#
# Copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Oct, 2011
# 
# R functions to do QCL mapping
#

#-- QCLscan.network interface --#
QCLscan.network <- function(genotypes, phenotypes, pheno.col, qcl.threshold = 0.4, max.depth=10, verbose = TRUE){
  if(missing(genotypes)) stop("genotypes are missing")
  if(missing(phenotypes)) stop("phenotypes are missing")
  if(missing(pheno.col)) stop("pheno.col missing")
  depth <- 1
  qcl_done <- QCLscan(genotypes,phenotypes,pheno.col)
  qcl_todo <- colnames(QCLprofiles(qcl_done,qcl.threshold,against="phenotypes"))
  while(depth <= max.depth){
    names_done <- unlist(lapply(qcl_done,attr,"name"))
    still_need_todo <- NULL
    for(x in qcl_todo){
      if(!(x %in% names_done)){
        still_need_todo <- c(still_need_todo,x)
      }
    }
    cat("Depth:",depth,"done:",length(names_done),"Todo:",length(qcl_todo),"/",length(still_need_todo),"\n")
    if(length(still_need_todo) > 0){
      ids <- which(colnames(phenotypes) %in% still_need_todo)
      tmp_scan <- QCLscan(genotypes,phenotypes,ids)
      qcl_todo <- colnames(QCLprofiles(tmp_scan,qcl.threshold,against="phenotypes"))
      qcl_done <- c(qcl_done,tmp_scan)
      depth <- depth+1
    }else{
      cat("Finished at depth:",depth,", scanned",length(qcl_done),"phenotypes\n")
      class(qcl_done) <- c(class(qcl_done),"QCLscan")
      return(qcl_done)
    }
  }
  cat("Reached max.depth, scanned",length(qcl_done),"phenotypes\n")
  class(qcl_done) <- c(class(qcl_done),"QCLscan")
  qcl_done
}

#-- R/qtl interface --#
QCLscan.network.cross <- function(cross, pheno.col, qcl.threshold = 0.4, max.depth=10, verbose = TRUE){
  if(missing(cross)) stop("cross is missing")
  if(missing(pheno.col)) stop("pheno.col missing")
  if(get(".has_rqtl", envir = .QclEnv)){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    genotypes <- qtl::pull.geno(cross)  
    QCLscan.network(genotypes, phenotypes, pheno.col=pheno.col, 
                    qcl.threshold=qcl.threshold,max.depth=max.depth,verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}
