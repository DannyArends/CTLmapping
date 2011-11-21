#
# qcl.profiles.R
#
# copyright (c) 2011 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Oct, 2010
# 
#

#Create the 2 possible QCL matrices
#Phenotypes versus markers
#Phenotypes versus phenotypes
QCLprofiles <- function(QCLscan, qcl.threshold=0.4, against = c("markers","phenotypes")){
  mymatrix <- NULL
  phenotypenames <- NULL
  if(against[1] == "markers"){
    for(p in 1:length(QCLscan)){
      mymatrix <- rbind(mymatrix,apply(abs(QCLscan[[p]]),2,function(x){length(which(x > qcl.threshold))}))
      phenotypenames <- c(phenotypenames,attr(QCLscan[[p]],"name"))
    }
    mymatrix <- matrix(unlist(mymatrix),length(QCLscan),ncol(QCLscan[[1]]))
    class(mymatrix) <- c(class(mymatrix),"P2Mmatrix")
    attr(mymatrix,"name") <- attr(QCLscan,"name")
    rownames(mymatrix) <- phenotypenames
    colnames(mymatrix) <- colnames(QCLscan[[1]])
    return(mymatrix)
  }
  if(against[1] == "phenotypes"){
    targets <- NULL
    for(p in 1:length(QCLscan)){
      targets <- unique(c(targets,unique(as.character(unlist(apply(abs(QCLscan[[p]]),2,function(x){names(which(x > qcl.threshold))}))))))
      phenotypenames <- c(phenotypenames,attr(QCLscan[[p]],"name"))
    }
    if(!is.null(targets)){
      mymatrix <- matrix(0,length(QCLscan),length(targets))
      colnames(mymatrix) <- targets
      for(p in 1:length(QCLscan)){
        current_table <- table(as.character(unlist(apply(abs(QCLscan[[p]]),2,function(x){names(which(x > qcl.threshold))}))))
        mymatrix[p,names(current_table)] <- current_table
      }
      class(mymatrix) <- c(class(mymatrix),"P2Pmatrix")
      attr(mymatrix,"name") <- attr(QCLscan,"name") 
      rownames(mymatrix) <- phenotypenames
      colnames(mymatrix) <- targets
      return(mymatrix)
    }
  }
}

QCLprofile <- function(QCL, qcl.threshold=0.4, verbose = TRUE){
  QCLscan <- list(QCL)
  class(QCLscan) <- c(class(QCLscan),"QCLscan")
  p2mm <- QCLprofiles(QCLscan, qcl.threshold, "markers")
  p2pm <- QCLprofiles(QCLscan, qcl.threshold, "phenotypes")
  sign_m <- colnames(p2mm)[which(p2mm > 0)]
  sign_p <- colnames(p2pm)[which(p2pm > 0)]
  if(verbose) cat("\t",substr(sign_m,1,6),"\n",sep="\t")
  result <- NULL
  for(phenoname in sign_p){
    myrow <- as.numeric(QCL[phenoname,sign_m]>qcl.threshold)
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
