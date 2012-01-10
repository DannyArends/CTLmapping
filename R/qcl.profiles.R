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
QCLprofiles <- function(QCLobject, significance = 0.05, against = c("markers","phenotypes"), verbose=FALSE){
  mymatrix <- NULL
  mynames <- NULL
  for(p in 1:length(QCLobject)){
    if(verbose) cat("Processing:",p,"from QCL to LOD\n")
    lod <- QCLtoLODvector(QCLobject[[p]], against)
    cat(length(lod),"\n")
    if(max(lod) > -log10(significance)){
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

QCLprofile <- function(QCLobject, pheno.col=1, significance = 0.05, verbose = TRUE){
  p2mm <- QCLprofiles(QCLobject, significance, "markers", verbose)
  p2pm <- QCLprofiles(QCLobject, significance, "phenotypes", verbose)
  sign_m <- colnames(p2mm)[which(apply(p2mm,2,function(x){any(x > -log10(significance))}))]
  sign_p <- colnames(p2pm)[which(apply(p2pm,2,function(x){any(x > -log10(significance))}))]
  if(verbose) cat("\t",substr(sign_m,1,6),"\n",sep="\t")
  result <- NULL
  for(phenoname in sign_p){
    myrow <- as.numeric(QCLobject[[pheno.col]]$l[phenoname,sign_m] > -log10(significance))
    if(verbose) cat(substr(phenoname,1,6),sum(myrow),myrow,"\n",sep="\t")
    result <- rbind(result, myrow)
  }
  if(!is.null(result)){
    if(verbose) cat("\t",apply(result,2,sum),"\n",sep="\t")
    result <- matrix(result,length(sign_p),length(sign_m))
    rownames(result) <- sign_p
    colnames(result) <- sign_m
    class(result) <- c(class(result),"QCLmatrix")
    return(result)
  }
}

# end of qcl.profiles.R
