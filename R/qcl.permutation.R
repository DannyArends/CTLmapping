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

#-- QCLpermute main function --#
QCLpermuteMC <- function(genotypes, phenotypes, pheno.col, n.perm=10, n.cores=2, directory="permutations", verbose=FALSE, ...){
  if(!has_snow()) stop("SNOW is not installed (or not yet loaded)")
  if(!file.exists(directory)) dir.create(directory)
  if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
  QCLperm <- vector("list",length(pheno.col))
  idx <- 1
  cl <- makeCluster(rep("localhost",n.cores))
  for(p in pheno.col){
    ss <- proc.time()
    cat("Starting multi core permutation for phenotype",p,"\n")
    parLapply(cl,as.list(1:n.perm),get("QCLpermute.internal"), genotypes=genotypes, phenotypes=phenotypes, pheno.col=p,directory=directory,verbose=verbose)
    QCLperm[[idx]] <- read.QCLpermute(directory, p, n.perm, verbose)
    idx <- idx+1
    el <- proc.time()
    cat("-",n.perm,"permutations took:",as.numeric(el[3]-ss[3]),"seconds.\n")
  }
  stopCluster(cl)
  invisible(QCLperm)
}

QCLpermute.internal <- function(perm, genotypes, phenotypes, pheno.col, directory="permutations", verbose=FALSE, ...){
  require(qcl)
  sl <- proc.time()
  if(verbose) cat("- Starting permutation",perm,"\n")
  phenotypes <- phenotypes[sample(nrow(phenotypes)),]
  write.table(QCLscan(genotypes, phenotypes, pheno.col, ...),file=paste(directory,"/Permutation_",pheno.col,"_",perm,".txt",sep=""))
  el <- proc.time()
  if(verbose) cat("- permutation",perm,"took:",as.numeric(el[3]-sl[3]),"seconds.\n")
}

#-- QCLpermute main function --#
QCLpermute <- function(genotypes, phenotypes, pheno.col, n.perm=10, n.cores=2, directory="permutations", verbose=FALSE, ...){
  if(has_snow()){
    require(snow)
    cat("- SNOW found using",n.cores,"processors for calculation\n")
    QCLpermuteMC(genotypes, phenotypes, pheno.col, n.perm=n.perm, n.cores=n.cores, directory=directory, verbose=verbose, ...)
  }else{
    if(!file.exists(directory)) dir.create(directory)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    QCLperm <- vector("list",length(pheno.col))
    idx <- 1
    for(p in pheno.col){
      ss <- proc.time()
      cat("Starting permutation for phenotype",p,"\n")
      for(x in 1:n.perm){
        QCLpermute.internal(x,genotypes,phenotypes,p,directory,verbose)
      }
      cat("-",n.perm,"permutations took:",as.numeric(el[3]-ss[3]),"seconds.\n")
      QCLperm[[idx]] <- read.QCLpermute(directory, p, n.perm, verbose)
      idx <- idx+1
    }
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

significance.QCLHotSpot <- function(QCLpermute,significant=212){
  significant <- print.QCLpermute(QCLpermute)[3]
  maximums <- lapply(QCLpermute,function(x){max(apply(x,1,function(x){sum(x > significant)}))})
  sorted <- sort(unlist(maximums))
  l <- length(sorted)
  values <- NULL
  valnames <- NULL
  for(x in c(.95,.99,.999)){
    values <- c(values,sorted[l*x])
    valnames <- c(valnames,paste((1-x)*100,"%"))
    cat((1-x)*100,"%\t",sorted[l*x],"\n")
  }
  names(values) <- valnames
  invisible(values)
}

#-- R/qtl interface --#
QCLpermute.cross <- function(cross, pheno.col, n.perm=10, n.cores=2, directory="permutations", verbose=FALSE, ...){
  if(missing(cross)) stop("cross is missing")
  if(missing(pheno.col)) stop("pheno.col missing")
  if(has_rqtl()){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    genotypes <- qtl::pull.geno(cross)
    QCLpermute(genotypes, phenotypes, pheno.col=pheno.col, n.perm=n.perm, n.cores=n.cores, directory=directory, verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}
