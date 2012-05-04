#
# ctl.utils.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Feb, 2012
# first written nov, 2010
# 
# check.genotypes, getRVM, lodscorestoscanone, getCorrelatedPhenotypes, gcLoop

check.genotypes <- function(genotypes, geno.enc=c(1,2), verbose=FALSE){
  if(verbose) cat(" - Genotypes: individuals=",nrow(genotypes),", markers=",ncol(genotypes),"\n",sep="")
  if(length(geno.enc)!=2) stop("argument 'geno.enc' length is incorrect, provide two genotype.values")
  
  toremove <- NULL
  idx <- 1
  checks <- apply(genotypes,2,function(geno){
    if(length(which(geno==geno.enc[1])) == 0 | length(which(geno==geno.enc[2])) == 0){
      if(verbose) cat("Severe: Empty group, removing marker",idx,"\n")
      toremove <<- c(idx,toremove)
    }
    idx <<- idx+1
  })
  if(length(toremove) > 0) cat(" - Genotypes: removing",length(toremove),"markers\n")
  toremove
}

#Create a matrix with row length = n.perms, filled with random numbers 1..n.rows
getRVM <- function(n.perms, n.rows){
  rvm <- NULL
  for(x in 1:n.perms){
    rvm <- rbind(rvm,sample(n.rows))
  }
  rvm
}

#Change any list of lodscores into a scanone object (only pre-req: length(lodscores)==sum(nmar(cross))
lodscorestoscanone <- function(cross,lodscores,traitnames = NULL){
  if(has_rqtl()){
    require(qtl)
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
    qtlprofile
  }else{
    warning(.has_rqtl_warnmsg)
  }
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
