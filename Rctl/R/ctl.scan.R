#
# ctl.scan.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Feb, 2012
# first written Jan, 2011
# 
# R functions to do CTL mapping
#

#-- CTLscan main function --#
CTLscan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), n.perms=100, conditions = NULL, geno.enc=c(1,2), verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  st <- proc.time()
  cat("Stage 0.0: Checking data\n")
  toremove <- check.genotypes(genotypes, geno.enc, verbose)
  genotypes[genotypes==geno.enc[1]] <- 0
  genotypes[genotypes==geno.enc[2]] <- 1
  genotypes[is.na(genotypes)]       <- -999
  phenotypes[is.na(phenotypes)]     <- -999
  results <- vector("list",length(pheno.col))
  stage <- 1
  cat("[INFO] ",nrow(phenotypes)," individuals, ",ncol(genotypes)," markers\n")

  if(!is.null(toremove)){
    cat("Stage 0.",stage,": Cleaning genotype data for mapping\n",sep="")
    warning(paste("Removing",length(toremove),"/",ncol(genotypes),"genotype markers")) 
    genotypes <- genotypes[,-toremove]
  }
  idx <- 1
  for(phe in pheno.col){
    results[[idx]] <- cCTLmap(genotypes, phenotypes, pheno.col=phe, n.perms = n.perms)
    class(results[[idx]]) <- c(class(results[[idx]]),"CTLscan")
    idx <- idx + 1
  }
  cat("Done after: ",(proc.time()-st)[3]," seconds\n")
  class(results) <- c(class(results),"CTLobject")
  results
}

cCTLmap <- function(genotypes, phenotypes, pheno.col=1, n.perms=100, verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  n.ind = nrow(genotypes); n.mar = ncol(genotypes); n.phe = ncol(phenotypes)

  ss <- proc.time()
	result <- .C("R_mapctl",as.integer(n.ind), as.integer(n.mar), as.integer(n.phe),
                    			as.integer(unlist(genotypes)), as.double(unlist(phenotypes)),
                          as.integer((pheno.col-1)), as.integer(n.perms),
                          dcor=as.double(rep(0,n.mar*n.phe)),
                          ctl=as.double(rep(0,n.mar*n.phe)), 
                          PACKAGE="ctl")
  res <- list()
  res$dcor <- matrix(result$dcor, n.mar, n.phe)
  res$ctl  <- matrix(result$ctl, n.mar, n.phe)
  rownames(res$dcor) <- colnames(genotypes); colnames(res$dcor) <- colnames(phenotypes)
  rownames(res$ctl)  <- colnames(genotypes); colnames(res$ctl)  <- colnames(phenotypes)
  ee <- proc.time()
  if(verbose) cat("  - CTLscan of",colnames(phenotypes)[pheno.col],"took",as.numeric(ee[3]-ss[3]),"seconds\n")
  invisible(res)
}

#-- R/qtl interface --#
CTLscan.cross <- function(cross, pheno.col, method = c("pearson", "pearsonordered", "kendall", "spearman"), n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(has_rqtl()){
    require(qtl)

    geno.enc <- NULL
    if(any(class(cross)=="bc") || any(class(cross)=="riself") || any(class(cross)=="risib")) geno.enc <- c(1,2)
    if(is.null(geno.enc)) stop("Class of the cross needs to be either: riself,risib or bc")
    rqtl_pheno <- qtl::pull.pheno(cross)
    rqtl_c <- which(colnames(rqtl_pheno) %in% c("sex","pgm")) #R/qtl adds additional columns sex and pgm
    if(length(rqtl_c) > 0) rqtl_pheno <- rqtl_pheno[,-rqtl_c]
    phenotypes <- apply(rqtl_pheno,2,as.numeric)         #R/qtl phenotypes are a data.frame (Need matrix)
    if(missing(pheno.col)) pheno.col <- 1:ncol(phenotypes)
    genotypes <- qtl::pull.geno(cross)
    CTLscan(genotypes = genotypes, phenotypes =  phenotypes, geno.enc=geno.enc, 
            pheno.col = pheno.col, method=method, n.perm=n.perm, n.cores=n.cores, 
            directory=directory, saveFiles=saveFiles, verbose=verbose)
  }else{
    warning(.has_rqtl_warnmsg)
  }
}

# end of ctl.scan.R
