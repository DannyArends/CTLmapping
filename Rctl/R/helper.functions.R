#
# helper.functions.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written May, 2012
# 
# Basic QC routines used in the examples of CTL analysis
# - Get a CTLobject's name
# - Remove the diagonal from a matrix
# - Color range for plots (Red-Black-Blue)
# - Chromosome edge locations from mapfile
# - Get the top-correlated metabolites

ctl.version    <- function(){ return(c(1,0,0)) }

ctl.names      <- function(CTLobject){ unlist(lapply(CTLobject, ctl.name)) }
ctl.qtlmatrix  <- function(CTLobject){ return(attr(CTLobject,"qtl")); }

ctl.name       <- function(CTLscan){ return(attr(CTLscan,"name")); }
ctl.dcormatrix <- function(CTLscan){ return(CTLscan$dcor); }
ctl.qtlprofile <- function(CTLscan){ return(CTLscan$qtl); }
ctl.ctlmatrix  <- function(CTLscan){ return(CTLscan$ctl); }

remove.diag    <- function(x){ return(x*lower.tri(x) + x*upper.tri(x)); }
redblue        <- function(){c(rgb(abs(seq(-2,-0,0.1))/2,0,0), rgb(0,0,seq(0.1,2,0.1)/2))}
whiteblack     <- function(){c("white",gray.colors(10)[10:1])}

dcor <- function(genotypes, phenotypes, marker=1, pheno1=1, pheno2=1, geno.enc=c(1,2), verbose = FALSE){
  idx1 <- which(genotypes[,marker] == geno.enc[1])
  idx2 <- which(genotypes[,marker] == geno.enc[2])
  c1 <- cor(phenotypes[idx1,pheno1], phenotypes[idx1,pheno2])
  c2 <- cor(phenotypes[idx2,pheno1], phenotypes[idx2,pheno2])
  dcor <- (c1-c2)^2
  if(verbose){
    cat("COR_1: ",c1,", COR_2: ",c2,"\n",sep="")
    cat("DCOR: ",dcor,"\n",sep="")
  }
  invisible(return(c(c1,c2,dcor)))
}

get.chr.edges <- function(mapinfo){
  unlist(lapply(unique(mapinfo[,1]),function(x){
    max(which(mapinfo[,1]==x));
  }))
}

top.correlated <- function(x){
  ret <- t(apply(remove.diag(x), 1, function(r) {
    id <- which.max(abs(r))
    return(c(names(r)[id], id, r[id]))
  }))
  colnames(ret) <- c("top.correlated", "id", "correlation")
  return(ret)
}

# end of helper.functions.R
