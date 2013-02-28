#
# ctl.utils.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written nov, 2010
# 
# check.genotypes, getRVM, lodscorestoscanone, getCorrelatedPhenotypes, gcLoop

check.genotypes <- function(genotypes, geno.enc=c(1,2), minAmount = 20, verbose=FALSE){
  if(verbose) cat("check.genotypes: ", ncol(genotypes)," markers, ", nrow(genotypes)," individuals\n")
  if(length(geno.enc) < 2) stop("argument 'geno.enc' length is incorrect, provide at least two genotype.values")
  
  toremove <- NULL
  idx <- 1
  geno.encNaN <- c(geno.enc, NaN, NA)
  checks <- apply(genotypes,2 , function(geno){
    #We need at least 3 markers of a certain genotype
    for(x in geno.enc){
      if(length(which(geno==x)) <= minAmount){
        if(verbose) cat("Severe: Small/Empty group", x ," (size:",length(which(geno==x)),"), removing marker",idx,"\n")
        toremove <<- c(toremove, idx)
      }
    }
    if(any((geno %in% geno.encNaN) == FALSE)){
      if(verbose) cat("Severe: Unknown genotypes, removing marker",idx,"\n")
      toremove <<- c(toremove, idx)
    }
    idx <<- idx+1
  })
  toremove <- unique(toremove)
  if(length(toremove) > 0){
    cat("check.genotypes: Removing", length(toremove),"/",ncol(genotypes), "markers\n")
    cat(toremove, "\n")
  }
  invisible(toremove)
}

#Create a matrix with row length = n.perms, filled with random numbers 1..n.rows
getRVM <- function(n.perms, n.rows){
  rvm <- NULL
  for(x in 1:n.perms){ rvm <- rbind(rvm,sample(n.rows)); }
  rvm
}

#Change any list of lodscores into a scanone object (only pre-req: length(lodscores)==sum(nmar(cross))
lodscorestoscanone <- function(cross,lodscores,traitnames = NULL){
  mymap <- qtl::pull.map(cross)
  n <- unlist(lapply(FUN=names,mymap))
  chr <- NULL
  if(!is.null(ncol(mymap[[1]]))){
    d <- as.numeric(unlist(lapply(mymap,FUN=function(x) {x[1,]})))
    for(i in 1:qtl::nchr(cross)){
      chr <- c(chr,rep(names(cross$geno)[i], ncol(mymap[[i]])))
    }
  }else{
    d <- as.numeric(unlist(mymap))
    for(i in 1:qtl::nchr(cross)){
      chr <- c(chr,rep(names(cross$geno)[i], length(mymap[[i]])))
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
  invisible(qtlprofile)
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

# end of ctl.utils.R
