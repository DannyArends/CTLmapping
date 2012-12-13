#
# ctl.lineplot.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends and Ritsert C. Jansen
# last modified Dec, 2012
# first written Dec, 2012
# 
# Line plot routines for CTL analysis
#

ctl.lineplot <- function(CTLobject, mapinfo, pheno.col, significance = 0.05, gap = 50, cex = 1, verbose = FALSE){
  if(missing(CTLobject) || is.null(CTLobject)) stop("argument 'CTLobject' is missing, with no default")
  if(missing(pheno.col)) pheno.col <- 1:length(CTLobject) 
  CTLobject  <- CTLobject[pheno.col]
  n.markers  <- nrow(CTLobject[[1]]$ctl)
  ctls       <- CTLnetwork(CTLobject, mapinfo, significance, verbose = verbose)
  if(is.null(ctls)) return() # No ctls found

  markerlocs <- cbind(seq(1,n.markers),rep(0,n.markers))
  total.l    <- n.markers  
  if(!missing(mapinfo)){
    markerlocs <- mapinfotomarkerlocs(mapinfo, gap, "line")
    total.l    <- max(markerlocs[,1])
  }
  xdim <- c(min(markerlocs[,1]),max(markerlocs[,1]))
  fromlocs <- total.l / (length(nfrom(ctls))+1)
  tolocs   <- total.l / (length(nto(ctls))+1)

  plot(xdim, c(-1.1, 1.1), t='n', axes = FALSE, xlab = "", ylab = "")
  points(markerlocs, pch=20, cex=(cex/2))
  for(x in 1:nrow(ctls)){ # Plot the ctls
    from <- c(which(nfrom(ctls) %in% ctls[x,1]) * fromlocs,  0.6)
    to   <- c(which(nto(ctls) %in% ctls[x,3]) * tolocs, -0.6)
    via  <- c(markerlocs[ctls[x,2]],  0.0)
    draw.spline(from, to, via, lwd=(ctls[x,4]/5)+1,col=ctls[x,1])
  } # All done now plot the trait elements
  for(x in 1:length(nfrom(ctls))){
    px <- which(nfrom(ctls) %in% nfrom(ctls)[x]) * fromlocs
    draw.element(px, 0.6, pheno.col[nfrom(ctls)[x]], cex=cex)
  }
  for(x in 1:length(nto(ctls))){
    px <- which(nto(ctls) %in% nto(ctls)[x]) * tolocs
    draw.element(px, -0.6, nto(ctls)[x], cex=cex)
  }
  box()
}

