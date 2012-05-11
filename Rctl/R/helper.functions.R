#
# helper.functions.R
#
# copyright (c) 2010-2011 Danny Arends and Ritsert C. Jansen
# last modified May, 2012
# first written May, 2012
# 
# Basic QC routines used in the examples of CTL analysis
# - Get a CTLobject's name
# - Remove the diagonal from a matrix
# - Color range for plots (Red-Black-Blue)
# - Chromosome edge locations from mapfile
# - Get the top-correlated metabolites

ctl.names      <- function(CTLobject){ unlist(lapply(CTLobject,function(x){return(attr(x$ctl,"name"));})) }
ctl.qtlmatrix  <- function(CTLobject){ return(attr(CTLobject,"qtl")); }

ctl.name       <- function(CTLscan){ return(attr(CTLscan$ctl,"name")); }
ctl.ctlmatrix  <- function(CTLscan){ return(CTLscan$ctl); }
ctl.qtlprofile <- function(CTLscan){ return(CTLscan$qtl); }
ctl.lodmatrix  <- function(CTLscan){ return(CTLscan$l); }

remove.diag <- function(x){ return(x*lower.tri(x) + x*upper.tri(x)); }
up <- function(){abs(seq(-2,-0,0.1))/2} 
dw <- function(){seq(0.1,2,0.1)/2}
redblue <- function(){c(rgb(up,0,0), rgb(0,0,dw))}
whiteblack <- function(){c("white",gray.colors(10)[10:1])}

get.chr.edges <- function(mapinfo){
  unlist(lapply(unique(mapinfo[,1]),function(x){max(which(mapinfo[,1]==x));}))
}

top.correlated <- function(x){
  ret <- t(apply(remove.diag(x),1,function(r){
    id <- which.max(abs(r))
    return(c(names(r)[id],id,r[id]))
  }))
  colnames(ret) <- c("top.correlated","id","correlation")
  return(ret)
}

