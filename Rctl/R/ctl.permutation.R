#
# ctl.permutations.R
#
# copyright (c) 2010 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Jan, 2012
# first written Nov, 2010
# 
# R functions to do permutation on CTL mapping
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#

#-- CTLpermuteMC function --#
CTLpermuteMC <- function(genotypes, phenotypes, geno.enc=c(1,2), pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(!has_snow()) stop("SNOW is not installed (or not yet loaded)")
  if(!file.exists(directory)) dir.create(directory)
  if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
  CTLperm <- vector("list",length(pheno.col))
  idx <- 1
  cl <- snow::makeCluster(rep("localhost",n.cores))
  for(p in pheno.col){
    ss <- proc.time()
    if(verbose) cat("  - Starting multi core permutation for phenotype",p,"\n")
    rvm <- getRVM(n.perm,nrow(genotypes))
    CTLperm[[idx]] <- unlist(snow::parLapply(cl,as.list(1:n.perm),get("CTLpermute.internal"), genotypes, phenotypes, geno.enc, p, method, rvm, directory, saveFiles, verbose,...))
    idx <- idx+1
    el <- proc.time()
    if(verbose) cat("  -",n.perm,"permutations took:",as.numeric(el[3]-ss[3]),"seconds.\n")
  }
  class(CTLperm) <- c(class(CTLperm),"CTLpermute")
  snow::stopCluster(cl)
  invisible(CTLperm)
}

CTLpermute.internal <- function(perm, genotypes, phenotypes, geno.enc=c(1,2), pheno.col, method = c("pearson", "kendall", "spearman"), rvm, directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  require(ctl)
  sl <- proc.time()
  if(verbose) cat("  - Starting permutation",perm,"\n")
  genotypes <- genotypes[rvm[perm,],]
  myperm <- CTLmapping(genotypes, phenotypes, geno.enc, pheno.col, method,verbose=verbose)
  if(saveFiles) write.table(myperm, file=paste(directory,"/Permutation_",pheno.col,"_",perm,".txt",sep=""))
  el <- proc.time()
  if(verbose) cat("  - Permutation",perm,"took:",as.numeric(el[3]-sl[3]),"seconds.\n")
  as.numeric(apply(abs(myperm),2,max))
}

#-- CTLpermute main function --#
CTLpermute <- function(genotypes, phenotypes, geno.enc=c(1,2), pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=10, n.cores=2, directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(has_snow() && n.cores > 1){
    require(snow)
    if(verbose) cat("  - SNOW found using",n.cores,"processors for calculation\n")
    CTLpermuteMC(genotypes, phenotypes, geno.enc, pheno.col, method, n.perm, n.cores, directory, saveFiles, verbose)
  }else{
    if(!file.exists(directory)) dir.create(directory)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    CTLperm <- vector("list",length(pheno.col))
    idx <- 1
    for(p in pheno.col){
      ss <- proc.time()
      if(verbose) cat("  - Starting permutation for phenotype",p,"\n")
      #Generate random numbers from a single thread, so we don't run into concurrency issues
      rvm <- getRVM(n.perm,nrow(genotypes))
      for(x in 1:n.perm){
        CTLperm[[idx]] <- c(CTLperm[[idx]],CTLpermute.internal(x, genotypes, phenotypes, geno.enc,p,method,rvm,directory,saveFiles,verbose))
      }
      el <- proc.time()
      if(verbose) cat("  -",n.perm,"permutations took:",as.numeric(el[3]-ss[3]),"seconds.\n")
      idx <- idx+1
    }
    class(CTLperm) <- c(class(CTLperm),"CTLpermute")
    invisible(CTLperm)
  }
}

#-- R/qtl interface --#
CTLpermute.cross <- function(cross, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=10, n.cores=2, geno.enc=c(1,2), directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(missing(cross)) stop("cross is missing")
  if(missing(pheno.col)) stop("pheno.col missing")
  if(has_rqtl()){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    genotypes <- qtl::pull.geno(cross)
    CTLpermute(genotypes, phenotypes, geno.enc=geno.enc, pheno.col=pheno.col, method=method, n.perm=n.perm, n.cores=n.cores, directory=directory, saveFiles = saveFiles, verbose=verbose,...)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}

# end of ctl.permutations.R
