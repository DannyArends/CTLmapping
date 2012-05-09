#
# read.dctl.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified May, 2012
# first written May, 2011
# 
# R functions to do transform CTL mapping scores to Pvalues and LOD
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#

read.dctl <- function(geno.file="genotypes.csv",pheno.file="phenotypes.csv",results="output/", verbose = TRUE){
  qtls.file <- paste(results,"qtls.txt",sep="")
  genotypes <- read.table(geno.file,row.names=1)
  phenotypes <- read.table(pheno.file,row.names=1)
  if(file.exists(qtls.file)){
    qtls <- read.table(qtls.file)
  }else{
    cat("[CTL] Output contains no QTL\n")
    genotypes <- genotypes[,(ncol(genotypes)-ncol(phenotypes)+1):ncol(genotypes)]
    qtls <- QTLscan(t(genotypes), t(phenotypes), verbose=verbose)$qtl  
  }
  cat("[CTL] Files loaded\n")
  rownames(qtls) <- rownames(phenotypes)
  colnames(qtls) <- rownames(genotypes)

  ctl_res <- list(nrow(phenotypes),mode="list")
  attr(ctl_res,"qtl") <- qtls

  for(idx in 1:nrow(phenotypes)){
    cat("[CTL] CTL scan ", idx," loaded\n",sep="")
    ctl <- read.table(paste(results,"ctl",(idx-1),".txt",sep=""))
    p <-   read.table(paste(results,"perms",(idx-1),".txt",sep=""))
    l <-   read.table(paste(results,"lodscores",(idx-1),".txt",sep=""))
    rownames(ctl) <- rownames(phenotypes)
    colnames(ctl) <- rownames(genotypes)
    rownames(l)   <- rownames(phenotypes)
    colnames(l)   <- rownames(genotypes)
    ctl_res[[idx]] <- list(4,mode="list")
    ctl_res[[idx]][[1]] <- ctl
    attr(ctl_res[[idx]][[1]],"name") <- rownames(phenotypes)[idx]
    class(ctl_res[[idx]][[1]]) <- c(class(ctl_res[[idx]][[1]]),"CTL")
    ctl_res[[idx]][[2]] <- qtls[idx,]
    ctl_res[[idx]][[3]] <- list()
    ctl_res[[idx]][[3]] <- as.numeric((unlist(p)))
    class(ctl_res[[idx]][[3]]) <- c(class(ctl_res[[idx]][[3]]),"CTLpermute")
    ctl_res[[idx]][[4]] <- l
    names(ctl_res[[idx]]) <- c("ctl","qtl","p","l")
    class(ctl_res[[idx]]) <- c(class(ctl_res),"CTLscan")
  }
  cat("[CTL] CTL object loaded\n")
  class(ctl_res) <- c(class(ctl_res),"CTLobject")
  return(ctl_res)
}
