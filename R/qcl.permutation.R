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
QCLpermuteMC <- function(genotypes, phenotypes, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=100, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(!has_snow()) stop("SNOW is not installed (or not yet loaded)")
  if(!file.exists(directory)) dir.create(directory)
  if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
  QCLperm <- vector("list",length(pheno.col))
  idx <- 1
  cl <- snow::makeCluster(rep("localhost",n.cores))
  for(p in pheno.col){
    ss <- proc.time()
    cat("  - Starting multi core permutation for phenotype",p,"\n")
    rvm <- getRVM(n.perm,nrow(genotypes))
    QCLperm[[idx]] <- unlist(snow::parLapply(cl,as.list(1:n.perm),get("QCLpermute.internal"), genotypes, phenotypes, p, method, rvm, genotype.values, directory, saveFiles, verbose,...))
    idx <- idx+1
    el <- proc.time()
    cat("  -",n.perm,"permutations took:",as.numeric(el[3]-ss[3]),"seconds.\n")
  }
  class(QCLperm) <- c(class(QCLperm),"QCLpermute")
  snow::stopCluster(cl)
  invisible(QCLperm)
}

QCLpermute.internal <- function(perm, genotypes, phenotypes, pheno.col, method = c("pearson", "kendall", "spearman"), rvm, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  require(qcl)
  sl <- proc.time()
  if(verbose) cat("  - Starting permutation",perm,"\n")
  genotypes <- genotypes[rvm[perm,],]
  myperm <- QCLmapping(genotypes, phenotypes, pheno.col, method, genotype.values,verbose=verbose)
  if(saveFiles) write.table(myperm, file=paste(directory,"/Permutation_",pheno.col,"_",perm,".txt",sep=""))
  el <- proc.time()
  if(verbose) cat("  - Permutation",perm,"took:",as.numeric(el[3]-sl[3]),"seconds.\n")
  as.numeric(max(abs(myperm)))
}

#-- QCLpermute main function --#
QCLpermute <- function(genotypes, phenotypes, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=10, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(has_snow() && n.cores > 1){
    require(snow)
    cat("  - SNOW found using",n.cores,"processors for calculation\n")
    QCLpermuteMC(genotypes, phenotypes, pheno.col, method, n.perm, n.cores, genotype.values, directory, saveFiles, verbose)
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
        QCLperm[[idx]] <- c(QCLperm[[idx]],QCLpermute.internal(x,genotypes,phenotypes,p,method,rvm,genotype.values,directory,saveFiles,verbose))
      }
      el <- proc.time()
      cat("  -",n.perm,"permutations took:",as.numeric(el[3]-ss[3]),"seconds.\n")
      idx <- idx+1
    }
    class(QCLperm) <- c(class(QCLperm),"QCLpermute")
    invisible(QCLperm)
  }
}

#-- R/qtl interface --#
QCLpermute.cross <- function(cross, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=10, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose=FALSE, ...){
  if(missing(cross)) stop("cross is missing")
  if(missing(pheno.col)) stop("pheno.col missing")
  if(has_rqtl()){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    genotypes <- qtl::pull.geno(cross)
    QCLpermute(genotypes, phenotypes, pheno.col=pheno.col, method=method, n.perm=n.perm, n.cores=n.cores, genotype.values=genotype.values, directory=directory, saveFiles = saveFiles, verbose=verbose,...)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}
