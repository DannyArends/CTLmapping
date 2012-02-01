#
# qcl.profiles.R
#
# copyright (c) 2011 Danny Arends and Ritsert C. Jansen
# last modified Jan, 2012
# first written Oct, 2010
# 
#

#Create the 2 possible QCL matrices
#Phenotypes versus markers
#Phenotypes versus phenotypes
QCLprofiles <- function(QCLobject, against = c("markers","phenotypes"), significance=0.05, verbose=FALSE, warn = TRUE){
  mymatrix <- NULL
  mynames <- NULL
  warn <- warn
  notice <- TRUE
  for(p in 1:length(QCLobject)){
    lod <- QCLtoLODvector(QCLobject[[p]], against)
    threshold <- -log10(significance)
    if(!is.nan(getPermuteThresholds(QCLobject[[p]])[1])){
      if(notice) cat("  - [Notice] Permutation available parameter 'significance' is ignored\n")
      threshold <- getPermuteThresholds(QCLobject[[p]])[1]
      notice <- FALSE
    }else{
      if(warn) cat("  - [Warning] Too few permutations, unable to find significant\n")
      warn <- FALSE
    }
    if(max(lod) > threshold){
      mymatrix <- rbind(mymatrix,lod)
      mynames <- c(mynames,attr(QCLobject[[p]]$qcl,"name"))  
    }
  }
  rownames(mymatrix) <- mynames
  if(against[1] == "phenotypes"){
    class(mymatrix) <- c(class(mymatrix),"P2Pmatrix")
  }else{
    class(mymatrix) <- c(class(mymatrix),"P2Mmatrix")
  }
  mymatrix
}

QCLprofile <- function(QCLobject, pheno.col=1, significance=0.05, verbose = TRUE){
  p2mm <- QCLprofiles(QCLobject, "markers", significance, verbose, FALSE)
  p2pm <- QCLprofiles(QCLobject, "phenotypes", significance, verbose, FALSE)
  threshold <- -log10(significance)
  if(!is.nan(getPermuteThresholds(QCLobject[[pheno.col]])[1])){
    cat("  - [Notice] Permutation available parameter 'significance' is ignored\n")
    threshold <- getPermuteThresholds(QCLobject[[pheno.col]])[1]
  }else{
    cat("  - [Warning] Too few permutations, using the significe parameter\n")
  }
  sign_m <- colnames(p2mm)[which(p2mm[pheno.col,] > threshold)]
  sign_p <- colnames(p2pm)[which(p2pm[pheno.col,] > threshold)]
  result <- NULL
  for(phenoname in sign_p){
    myrow <- as.numeric(QCLobject[[pheno.col]]$l[phenoname,sign_m] > threshold)
    result <- rbind(result, myrow)
  }
  if(!is.null(result)){
    result <- matrix(result,length(sign_p),length(sign_m))
    rownames(result) <- sign_p
    colnames(result) <- sign_m
    class(result) <- c(class(result),"QCLmatrix")
    result <- result[,-which(apply(result,2,sum)==0)]
    if(verbose){
      cat("\t",substr(colnames(result),1,6),"\n",sep="\t")
      for(x in 1:nrow(result)){
        cat(substr(rownames(result)[x],1,6),sum(result[x,]),result[x,],"\n",sep="\t")
      }
    }
    return(result)
  }
}

# end of qcl.profiles.R
