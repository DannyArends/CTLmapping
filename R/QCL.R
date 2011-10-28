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

QCLscan.network <- function(genotypes, phenotypes, pheno.col = 1, qcl.threshold = 0.3, m.depth=5, verbose = TRUE){
  todo <- pheno.col
  at_my_depth <- 1
  depth <- 1
  cnt <- 1
  results <- vector("list",1)
  while(depth <= m.depth){
    at_next_depth <- 0
    while(at_my_depth > 0){
      x <- todo[1]
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
      myprofile <- QCLprofiles(results,qcl.threshold=qcl.threshold,against="phenotypes")
      if(dim(myprofile)[2]!=0){
        additional <- which(colnames(phenotypes) %in% colnames(myprofile))
        done <- which(colnames(phenotypes)[additional] %in% unlist(lapply(results,attr,"name")))
        if(!is.na(done&&1))additional <- additional[-done]
        todo <- c(todo,additional)
        todo <- todo[-1]
        if(verbose){
          cat("At",depth,", added ",length(additional),"phenotypes\n")
        }
        at_next_depth <- at_next_depth + length(additional)
      }
      at_my_depth <- at_my_depth-1
      cnt <- cnt+1
      if(at_my_depth==0 && depth == m.depth){ 
      }else{
        length(results) <- length(results) +1
      }
    }
    at_my_depth <- at_next_depth
    depth <- depth +1
  }
  class(results) <- c(class(results),"QCLscan")
  results
}
