#
# qcl.utils.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Jun, 2011
# first written nov, 2010
# 
# Plotting routines for QCL analysis
#


QCLscanToProfile <- function(QCL, qcl.threshold=0.1){
  apply(QCL,2,function(x){length(which(x > qcl.threshold))})
}

#Change any list of lodscores into a scanone object (only pre-req: length(lodscores)==sum(nmar(cross))
lodscorestoscanone <- function(cross,lodscores,traitnames = NULL){
  n <- unlist(lapply(FUN=names,pull.map(cross)))
  chr <- NULL
  if(!is.null(ncol(pull.map(cross)[[1]]))){
    d <- as.numeric(unlist(lapply(pull.map(cross),FUN=function(x) {x[1,]})))
    for(i in 1:nchr(cross)){
      chr <- c(chr,rep(names(cross$geno)[i], ncol(pull.map(cross)[[i]])))
    }
  }else{
    d <- as.numeric(unlist(pull.map(cross)))
    for(i in 1:nchr(cross)){
      chr <- c(chr,rep(names(cross$geno)[i], length(pull.map(cross)[[i]])))
    }
  }
  qtlprofile <- cbind(chr,d,lodscores)
  qtlprofile <- as.data.frame(qtlprofile)
  qtlprofile[,1] <- chr
  qtlprofile[,2] <- as.numeric(d)
  if(!is.null(ncol(lodscores))){
    for(x in 1:ncol(lodscores)){
      qtlprofile[,2+x] <- as.numeric(lodscores[,x])
    }
    traitnames = paste("lod",1:ncol(lodscores))
  }else{
     qtlprofile[,3] <- as.numeric(lodscores)
     traitnames = "lod"
  }
  rownames(qtlprofile) <- n
  colnames(qtlprofile) <- c("chr","cM",traitnames)
  class(qtlprofile) <- c("scanone", "data.frame")
  qtlprofile
}

getCorrelatedPhenotypes <- function(cross, pheno.col=1, top=200, min_threshold=0.4){
  new_cross <- cross
  phenotypes <- apply(pull.pheno(cross),2,as.numeric)
  target <- phenotypes[,pheno.col]
  cors <- apply(phenotypes,2,function(x){cor(x,target,use="pair")})
  threshold <- min_threshold
  significant_cors <- which(cors >= threshold)
  while(length(significant_cors) > top){
    threshold = threshold + 0.001
    significant_cors <- which(cors >= threshold)
  }
  new_cross$pheno <- as.data.frame(phenotypes[,significant_cors])
  cat("Selected:",length(significant_cors),"phenotypes with correlation >=",threshold,"\n")
  new_cross
}

gcLoop <- function(verbose = FALSE){
  p_usage <- gc()[2,3]
  n_usage <- gc()[2,3]
  while(n_usage < p_usage){
    p_usage = n_usage
    n_usage <- gc()[2,3]
    if(verbose) cat("GCloop ",n_usage," ",p_usage,"\n")
  }
}
