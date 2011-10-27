#
# qcl.significance.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified feb, 2011
# first written nov, 2010
# 
# R functions to do permutation on QCL mapping
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#

QCLpermute <- function(cross, pheno.col, n.perm=10, directory="permutations", verbose=TRUE, ...){
  if(!file.exists(directory)) dir.create(directory)
  if(missing(pheno.col)) pheno.col <- 1:nphe(cross)
  for(p in pheno.col){
    cat("Starting permutation for phenotype",p,"\n")
    for(x in 1:n.perm){
      sl <- proc.time()
      if(verbose) cat("- Starting permutation",x,"/",n.perm,"\n")
		  cross$pheno <- cross$pheno[sample(nind(cross)),]
      write.table(QCLscan(cross, pheno.col=p, ...),file=paste(directory,"/Permutation_",x,".txt",sep=""))
      el <- proc.time()
      if(verbose) cat("- Permutation",x,"took:",as.numeric(el[3]-sl[3]),"Seconds.\n")
    }
  }
  QCLpermute <- read.QCLpermute(directory,n.perm)
  invisible(QCLpermute)
}

read.QCLpermute <- function(directory="permutations", n.perm){
  files <- dir(directory)
  if(missing(n.perm)) n.perm <- length(files)
  if(n.perm < 1) stop(paste("No permutation files found in:",directory))
  QCLpermute <- vector("list", n.perm)
  for(x in 1:n.perm){
    cat("Trying to read:",paste(directory,"/Permutation_",x,".txt\n",sep=""))
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
