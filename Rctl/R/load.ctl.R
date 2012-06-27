load.ctl <- function(genotypes = "ngenotypes.txt", phenotypes = "nphenotypes.txt", output = "ctlout", n.pheno=250, verbose = FALSE){
  require(ctl)
  genotypes = "ngenotypes.txt"; phenotypes = "nphenotypes.txt"; output = "ctlout"; verbose = TRUE
  phenodata <- read.csv(phenotypes,sep="\t",header=FALSE)
  genodata  <- read.csv(genotypes,sep="\t",header=FALSE)
  singleperms <- FALSE
  results <- vector("list",nrow(phenodata))
  if(file.exists(paste(output,"/qtls.txt",sep=""))){
    if(verbose) cat("QTL results found\n")
    qtldata <- read.csv(paste(output,"/qtls.txt",sep=""),sep="\t",header=FALSE)
    rownames(qtldata) <- phenodata[,1]
    colnames(qtldata) <- genodata[,1]
    for(idx in 1:nrow(phenodata)){ results[[idx]]$qtl <- qtldata[idx, ] }
  }else{
    cat("Warning no QTL results found (-qtl)\n")
  }
  for(x in 0:min((nrow(phenodata)-1),n.pheno)){
    idx <- (x+1)
    if(file.exists(paste(output,"/ctl",x,".txt",sep=""))){
      if(verbose) cat("CTLs found for ",idx,"/",nrow(phenodata)," ",phenodata[idx,1],"\n")
      results[[idx]]$ctl <-  read.csv(paste(output,"/ctl",x,".txt",sep=""),sep="\t",header=FALSE)
      rownames(results[[idx]]$ctl) <- phenodata[,1]
      colnames(results[[idx]]$ctl) <- genodata[,1]
      class(results[[idx]]$ctl) <- c(class(results[[idx]]$ctl),"CTL")
      attr(results[[idx]]$ctl,"name") <-  phenodata[idx,1]
    }else{
      stop("No CTLS found for:", phenodata[idx,1],"\n")
    }
    if(!singleperms && file.exists(paste(output,"/perms",x,".txt",sep=""))){
      results[[idx]]$p <- apply(read.csv(paste(output,"/perms",x,".txt",sep=""),sep="\t",header=FALSE),1,max)
      class(results[[idx]]$p) <- c(class(results[[idx]]$p),"CTLpermute")
    }else{
      if(x==1) cat("Results are 'single permutation mode' (-sp) \n")
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

load.lude <- function(){
  genotypes = "ngenotypes.txt"
  phenotypes = "nphenotypes.txt"
  output = "ctlout"
  memory.limit(2000)
  load.ctl(genotypes, phenotypes, output)
}
