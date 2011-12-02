#
# qcl.utils.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Jun, 2011
# first written nov, 2010
# 
# nchr

getRVM <- function(n.perms, n.rows){
  rvm <- NULL
  for(x in 1:n.perms){
    rvm <- rbind(rvm,sample(n.rows))
  }
  rvm
}

#Change any list of lodscores into a scanone object (only pre-req: length(lodscores)==sum(nmar(cross))
lodscorestoscanone <- function(cross,lodscores,traitnames = NULL){
  if(get(".has_rqtl", envir = .QclEnv)){
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

getCorrelatedPhenotypes <- function(cross, pheno.col=1, top=200, min_threshold=0.4){
  if(get(".has_rqtl", envir = .QclEnv)){
    require(qtl)
    new_cross <- cross
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
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
