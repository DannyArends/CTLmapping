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
QCLpermute <- function(genotypes, phenotypes, pheno.col, n.perm=10, directory="permutations", verbose=FALSE, ...){
  if(!file.exists(directory)) dir.create(directory)
  if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
  QCLperm <- vector("list",length(pheno.col))
  idx <- 1
  for(p in pheno.col){
    ss <- proc.time()
    cat("Starting permutation for phenotype",p,"\n")
    for(x in 1:n.perm){
      sl <- proc.time()
      if(verbose) cat("- Starting permutation",x,"/",n.perm,"\n")
		  phenotypes <- phenotypes[sample(nrow(phenotypes)),]
      write.table(QCLscan(genotypes, phenotypes, pheno.col=p, ...),file=paste(directory,"/Permutation_",x,".txt",sep=""))
      el <- proc.time()
      if(verbose) cat("- Permutation",x,"took:",as.numeric(el[3]-sl[3]),"Seconds.\n")
    }
    cat("-",x,"Permutations took:",as.numeric(el[3]-ss[3]),"Seconds.\n")
    QCLperm[[idx]] <- read.QCLpermute(directory,n.perm,verbose)
    idx <- idx+1
  }
  invisible(QCLperm)
}

read.QCLpermute <- function(directory="permutations", n.perm, verbose){
  files <- dir(directory)
  if(missing(n.perm)) n.perm <- length(files)
  if(n.perm < 1) stop(paste("No permutation files found in:",directory))
  QCLpermute <- vector("list", n.perm)
  for(x in 1:n.perm){
    if(verbose) cat("Trying to read:",paste(directory,"/Permutation_",x,".txt\n",sep=""))
    QCLpermute[[x]]  <- read.table(paste(directory,"/Permutation_",x,".txt",sep=""))
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
QCLpermute.cross <- function(cross, pheno.col, n.perm=10, directory="permutations", verbose=TRUE, ...){
  if(missing(cross)) stop("cross is missing")
  if(missing(pheno.col)) stop("pheno.col missing")
  if(get(".has_rqtl", envir = .QclEnv)){
    require(qtl)
    phenotypes <- apply(qtl::pull.pheno(cross),2,as.numeric)
    genotypes <- qtl::pull.geno(cross)
    QCLpermute(genotypes, phenotypes,pheno.col=pheno.col,n.perm=n.perm,directory=directory,verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}
