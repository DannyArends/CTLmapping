load.ctl <- function(genotypes = "ngenotypes.txt", phenotypes = "nphenotypes.txt", output = "ctlout", from=1, to=250, verbose = FALSE){
  require(ctl)
  phenodata <- read.csv(phenotypes,sep="\t",header=FALSE)
  genodata  <- read.csv(genotypes,sep="\t",header=FALSE)
  singleperms <- FALSE
  results <- vector("list",(to-(from-1)))

  if(file.exists(paste(output,"/qtls.txt",sep=""))){
    if(verbose) cat("QTL results found\n")
    qtldata <- read.csv(paste(output,"/qtls.txt",sep=""),sep="\t",header=FALSE)
    rownames(qtldata) <- phenodata[,1]
    colnames(qtldata) <- genodata[,1]
    for(idx in from:min(nrow(phenodata),to)){ 
      cat(idx-(from-1),"\n")
      results[[idx-(from-1)]]$qtl <- qtldata[idx, ] 
    }
  }else{
    cat("Warning no QTL results found (-qtl)\n")
  }
  for(x in (from-1):min((nrow(phenodata)-1), (to-1))){
    idx <- (x+1)-(from-1) # D2.0 counts from 0, R from 1
    cat(idx," ", x," ",phenodata[(x+1),1],"\n")
    if(file.exists(paste(output,"/ctl",x,".txt",sep=""))){
      if(verbose) cat("CTLs found for ",idx,"/",min((nrow(phenodata)), to - from)," ",phenodata[idx,1],"\n")
      results[[idx]]$ctl <-  read.csv(paste(output,"/ctl",x,".txt",sep=""),sep="\t",header=FALSE)
      rownames(results[[idx]]$ctl) <- phenodata[,1]
      colnames(results[[idx]]$ctl) <- genodata[,1]
      class(results[[idx]]$ctl) <- c(class(results[[idx]]$ctl),"CTL")
      attr(results[[idx]]$ctl,"name") <-  phenodata[(x+1),1]
      
    }else{
      stop("No CTLS found for:", phenodata[x,1],"\n")
    }
    if(!singleperms && file.exists(paste(output,"/perms",x,".txt",sep=""))){
      results[[idx]]$p <- apply(read.csv(paste(output,"/perms",x,".txt",sep=""),sep="\t",header=FALSE),1,max)
      class(results[[idx]]$p) <- c(class(results[[idx]]$p),"CTLpermute")
    }else{
      if(x==1) cat("Results are 'single permutation mode' (-sp) \n")
      if(from > 1 && !singleperms) results[[1]]$p <- apply(read.csv(paste(output,"/perms0.txt",sep=""),sep="\t",header=FALSE),1,max)
      singleperms <- TRUE
      results[[idx]]$p <- results[[1]]$p
    }
    if(verbose) cat("Redo of LOD transformation using GPD distribution\n")
    results[[idx]]$l <- toLod(results[[idx]], FALSE, FALSE)
    class(results[[idx]]) <- c(class(results[[idx]]),"CTLscan")
  }
  class(results) <- c(class(results),"CTLobject")
  results
}

#load.lude <- function(){
  memory.limit(2000)
  setwd("e:/gbic/ludectl")
 
  resN3 <- load.ctl("n3genotypes.txt", "n3phenotypes.txt", "ctln3out",251,500)
   
  aa <- CTLprofiles(resN3,signi=0.01)
  write.table(aa,file="ctls_n3_500.txt",col.names=TRUE, append=FALSE, sep="\t")
#}
