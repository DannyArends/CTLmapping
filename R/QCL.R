#
# qcl.R
#
# copyright (c) 2010 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2011
# first written Oct, 2011
# 
# R functions to do QCL mapping
#

QCLscan.cross <- function(cross, pheno.col, verbose = FALSE){
  if(.has_rqtl){
    require(qtl)
    phenotypes <- apply(pull.pheno(cross),2,as.numeric)
    genotypes <- pull.geno(cross)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
  
    QCLscan(genotypes, phenotypes,pheno.col=pheno.col,verbose=verbose)
  }else{
    stop(.has_rqtl_warnmsg)
  }
}

QCLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), verbose = FALSE){
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

QCLscan.network <- function(genotypes, phenotypes, pheno.col = 1, qcl.threshold = 0.3, max.depth=10, verbose = TRUE){
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
