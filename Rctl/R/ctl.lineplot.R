#
# ctl.lineplot.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Dec, 2012
# first written Dec, 2012
# 
# Line plot routines for CTL analysis
#

ctl.line <- function(){
  n.mar <- 100
  plot(c(0, n.mar), c(-1.1, 1.1), t='n', axes = FALSE, xlab = "", ylab = "")
  points(cbind(seq(1,n.mar),rep(0,n.mar)), pch=20, cex=.1)
  nfrom <- function(ctls){ return(unique(ctls[,1])) }
  nto   <- function(ctls){ return(unique(ctls[,3])) }
  for(x in 1:nrow(ctls)){
    from <- c(which(nfrom(ctls) %in% ctls[x,1]) * (n.mar/length(nfrom(ctls))),  0.6)
    to   <- c(which(nto(ctls) %in% ctls[x,3]) * (n.mar/length(nto(ctls))), -0.6)
    via  <- c(ctls[x,2],  0.0)
    draw.spline(from, to, via, lwd=ctls[x,4]/5,col=ctls[x,1])
  }
  for(x in 1:length(nfrom(ctls))){
    px <- which(nfrom(ctls) %in% nfrom(ctls)[x]) * (n.mar/length(nfrom(ctls)))
    points(px, 0.6, cex=3, pch=20,col="white")
    points(px, 0.6, cex=3)
    text(px, 0.6, nfrom(ctls)[x]) 
  }
  for(x in 1:length(nto(ctls))){
    px <- which(nto(ctls) %in% nto(ctls)[x]) * (n.mar/length(nto(ctls)))
    points(px, -0.6, cex=3, pch=20,col="white")
    points(px, -0.6, cex=3)
    text(px, -0.6, nto(ctls)[x]) 
  }
}

