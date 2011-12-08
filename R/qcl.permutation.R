#
# qcl.permutations.R
#
# copyright (c) 2010 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2011
# first written nov, 2010
# 
# R functions to do permutation on QCL mapping
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#

#-- QCLpermuteMC function --#
QCLpermuteMC <- function(genotypes, phenotypes, pheno.col, n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(!has_snow()) stop("SNOW is not installed (or not yet loaded)")
  if(!file.exists(directory)) dir.create(directory)
  if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
  QCLperm <- vector("list",length(pheno.col))
  idx <- 1
  cl <- makeCluster(rep("localhost",n.cores))
  for(p in pheno.col){
    ss <- proc.time()
    cat("  - Starting multi core permutation for phenotype",p,"\n")
    rvm <- getRVM(n.perm,nrow(genotypes))
    QCLperm[[idx]] <- unlist(parLapply(cl,as.list(1:n.perm),get("QCLpermute.internal"), genotypes, phenotypes, p, rvm, directory, saveFiles, verbose))
    class(QCLperm[[idx]]) <- c(class(QCLperm[[idx]]),"QCLpermute")
    idx <- idx+1
    el <- proc.time()
    cat("  -",n.perm,"permutations took:",as.numeric(el[3]-ss[3]),"seconds.\n")
  }
  stopCluster(cl)
  invisible(QCLperm)
}

QCLpermute.internal <- function(perm, genotypes, phenotypes, pheno.col, rvm, directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  require(qcl)
  sl <- proc.time()
  if(verbose) cat("  - Starting permutation",perm,"\n")
  genotypes <- genotypes[rvm[perm,],]
  perm <- QCLmapping(genotypes, phenotypes, pheno.col, ...)
  if(saveFiles) write.table(perm, file=paste(directory,"/Permutation_",pheno.col,"_",perm,".txt",sep=""))
  el <- proc.time()
  if(verbose) cat("  - Permutation",perm,"took:",as.numeric(el[3]-sl[3]),"seconds.\n")
  as.numeric(max(abs(perm)))
}

#-- QCLpermute main function --#
QCLpermute <- function(genotypes, phenotypes, pheno.col, n.perm=10, n.cores=2, directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(has_snow() && n.cores > 1){
    require(snow)
    cat("  - SNOW found using",n.cores,"processors for calculation\n")
    QCLpermuteMC(genotypes, phenotypes, pheno.col, n.perm=n.perm, n.cores=n.cores, directory, saveFiles, verbose, ...)
  }else{
    if(!file.exists(directory)) dir.create(directory)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    QCLperm <- vector("list",length(pheno.col))
    idx <- 1
    for(p in pheno.col){
      ss <- proc.time()
      cat("  - Starting permutation for phenotype",p,"\n")
      #Generate random number from a single thread, so we don't run into concurrency issues
      rvm <- getRVM(n.perm,nrow(genotypes))
      for(x in 1:n.perm){
        QCLperm[[idx]] <- c(QCLperm[[idx]],QCLpermute.internal(x,genotypes,phenotypes,p,rvm,directory,saveFiles,verbose))
      }
      el <- proc.time()
      cat("  -",n.perm,"permutations took:",as.numeric(el[3]-ss[3]),"seconds.\n")
      idx <- idx+1
    }
    class(QCLperm) <- c(class(QCLperm),"QCLpermute")
    invisible(QCLperm)
  }
}

read.QCLpermute <- function(directory="permutations", pheno.col=1, n.perm, verbose){
  files <- dir(directory)
  if(missing(n.perm)) n.perm <- length(files)
  if(n.perm < 1) stop(paste("No permutation files found in:",directory))
  QCLpermute <- vector("list", n.perm)
  for(x in 1:n.perm){
    if(verbose) cat("Trying to read:",paste(directory,"/Permutation_",pheno.col,"_",x,".txt\n",sep=""))
    QCLpermute[[x]]  <- read.table(paste(directory,"/Permutation_",pheno.col,"_",x,".txt",sep=""))
  }
  class(QCLpermute) <- c(class(QCLpermute),"QCLpermute")
  invisible(QCLpermute)
}

#-- R/qtl interface --#
QCLpermute.cross <- function(cross, pheno.col, n.perm=10, n.cores=2, directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(missing(cross)) stop("cross is missing")
  if(missing(pheno.col)) stop("pheno.col missing")
  if(has_rqtl()){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    genotypes <- qtl::pull.geno(cross)
    QCLpermute(genotypes, phenotypes, pheno.col=pheno.col, n.perm=n.perm, n.cores=n.cores, directory=directory, saveFiles = saveFiles, verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}
