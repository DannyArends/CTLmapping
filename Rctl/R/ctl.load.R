#
# ctl.load.R
#
# copyright (c) 2012 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Jul, 2012
# first written Jun, 2012
# 
# R functions to do load CTL mapping donw with the D version
#

#-- ctl.load function --#
ctl.load <- function(genotypes = "ngenotypes.txt", phenotypes = "nphenotypes.txt", output = "ctlout", from=1, to, verbose = FALSE){
  require(ctl)
  phenodata   <- read.csv(phenotypes,sep="\t",header=FALSE)
  genodata    <- read.csv(genotypes,sep="\t",header=FALSE)
  singleperms <- FALSE
  if(missing(to)){
    to <- nrow(phenodata)
  }
  if(verbose) cat("Attempting to load ", to," CTL scans\n")
  results <- vector("list",(to-(from-1)))
  if(file.exists(paste(output,"/qtls.txt",sep=""))){
    if(verbose) cat("QTL results found\n")
    qtldata <- read.csv(paste(output,"/qtls.txt",sep=""),sep="\t",header=FALSE)
    rownames(qtldata) <- phenodata[,1]
    colnames(qtldata) <- genodata[,1]
    for(idx in from:min(nrow(phenodata),to)){ 
      results[[idx-(from-1)]]$qtl <- qtldata[idx, ] 
    }
  }else{
    cat("No QTL results found, use the --qtl option to also scan for QTLs\n")
  }
  for(x in (from-1):min((nrow(phenodata)-1), (to-1))){
    idx <- (x+1)-(from-1) # D2.0 counts from 0, R from 1
    #if(verbose) cat("Trying to load:", idx," ", x,":",phenodata[(x+1),1],"\n")
    if(file.exists(paste(output,"/ctl",x,".txt",sep=""))){
      if(verbose) cat("CTLs found for ",idx,"/",min((nrow(phenodata)), to - (from-1))," ",phenodata[idx,1],"\n")
      results[[idx]]$ctl <-  read.csv(paste(output,"/ctl",x,".txt",sep=""),sep="\t",header=FALSE)
      rownames(results[[idx]]$ctl)    <- phenodata[,1]
      colnames(results[[idx]]$ctl)    <- genodata[,1]
      class(results[[idx]]$ctl)       <- c(class(results[[idx]]$ctl),"CTL")
      attr(results[[idx]]$ctl,"name") <-  phenodata[(x+1),1]
    }else{
      stop("ERROR: Missing file No CTLs found for:", phenodata[x,1],"\n")
    }
    if(!singleperms && file.exists(paste(output,"/perms",x,".txt",sep=""))){
      results[[idx]]$perms        <- apply(read.csv(paste(output,"/perms",x,".txt",sep=""),sep="\t",header=FALSE),1,max)
      class(results[[idx]]$perms) <- c(class(results[[idx]]$perms),"CTLpermute")
    }else{
      if(x==1) cat("[Warning] results are from 'single permutation mode' (--sp) NOTE: This might induce many false positives \n")
      if(from > 1 && !singleperms) results[[1]]$perms <- apply(read.csv(paste(output,"/perms0.txt",sep=""),sep="\t",header=FALSE),1,max)
      singleperms       <- TRUE
      results[[idx]]$perms  <- results[[1]]$perms
    }
    if(verbose) cat("Recalculating LOD score using permutations, using GPD distribution for extreme scores\n")
    results[[idx]]$ctl    <- toLod(results[[idx]], FALSE, FALSE)
    class(results[[idx]]) <- c(class(results[[idx]]),"CTLscan")
  }
  class(results) <- c(class(results),"CTLobject")
  results
}
