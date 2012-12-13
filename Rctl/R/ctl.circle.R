#
# ctl.circle.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Dec, 2012
# first written Dec, 2012
# 
# Circle plot routines for CTL analysis
#

draw.spline <- function(cn1, cn2, via = c(0,0), lwd = 1, col="blue", ...){
  x <- c(cn1[1], via[1], cn2[1])
  y <- c(cn1[2], via[2], cn2[2])
  invisible(xspline(x, y, shape=0, lwd=lwd, border=col,...))
}

draw.element <- function(x, y, title, cex=3, bg.col = "white", border.col="black"){
  points(c(x, y), cex=cex, pch=20, col=bg.col)
  points(c(x, y), cex=cex, col= border.col)
  text(x, y, title)
}

circle.loc <- function(nt, size = 1.0){
  locs <- matrix(nrow = nt, ncol = 2)
  phi  <- seq(0, 2 * pi, length = (nt+1))
  complex.circle <- complex(modulus = 1, argument = phi)
  for(j in 1:nt){ locs[j, ] <- c(Im(complex.circle[j]), Re(complex.circle[j])); }
  invisible(locs * size)
}

ctls <- rbind(c(1, 4, 3, 3.2), c(1, 5, 3, 5.6), c(2, 6, 15, 4.3), c(2, 78, 7, 6.4), c(2, 30, 8, 15.0),
              c(1, 2, 7, 6.4), c(15, 21, 7, 15.0), c(15, 20, 2, 6.3), c(5, 40, 8, 6.4), c(5, 5, 9, 7.3))

nfrom <- function(ctls){ return(unique(ctls[,1])) }
nto   <- function(ctls){ return(unique(ctls[,3])) }

ctl.circle <- function(CTLobject, mapinfo, pheno.col, significance = 0.05, gap = 20, verbose = FALSE){
  if(missing(pheno.col)) pheno.col <- 1:length(CTLobject) 
  CTLobject  <- CTLobject[pheno.col]
  n.markers  <- nrow(CTLobject[[1]]$ctl)
  ctls       <- CTLnetwork(CTLobject, mapinfo, significance, verbose = verbose)
  if(is.null(ctls)) return()
  n.chr      <- 1
  markerlocs <- circle.loc(n.markers, 0.7)
  if(!missing(mapinfo)){
    n.chr       <- unique(mapinfo[,1])
    chr.lengths <- NULL
    for(chr in 1:length(n.chr)){
      ll <- mapinfo[lapply(unique(mapinfo[,1]),function(x){which(x==mapinfo[,1])})[[chr]],]
      chr.lengths <- c(chr.lengths, max(as.numeric(ll[,2]))-min(as.numeric(ll[,2])))
    }
    total.l <- ceiling(sum(chr.lengths) + (gap*length(n.chr)))
    cmmap   <- circle.loc(total.l, 0.7)
    markerlocs <- NULL
    for(x in 1:nrow(mapinfo)){
      m.chr <- which(unique(mapinfo[,1]) %in% mapinfo[x,1])
      m.loc <- ceiling(mapinfo[x,2])
      if(m.chr > 1){ m.loc <- m.loc + sum(chr.lengths[1:(m.chr-1)]) + ((m.chr-1)*gap); }
      markerlocs <- rbind(markerlocs, cmmap[ceiling(m.loc),])
    }
  }
  fromtlocs  <- circle.loc(length(nfrom(ctls)), 1.0)
  totlocs    <- circle.loc(length(nto(ctls)), 0.4)
  plot(c(-1.1, 1.1), c(-1.1, 1.1), type = "n", axes = FALSE, xlab = "", ylab = "")
  points(markerlocs, pch=20, cex=1)
  for(x in 1:nrow(ctls)){
    from <- fromtlocs[which(nfrom(ctls) %in% ctls[x,1]),]
    to   <- totlocs[which(nto(ctls) %in% ctls[x,3]),]
    via  <- markerlocs[ctls[x,2],]
    draw.spline(from, to, via, lwd=(ctls[x,4]/5)+1,col=ctls[x,1])
  }
  points(fromtlocs, cex=4, pch=20,col="white")
  points(fromtlocs, cex=4)
  for(x in 1:nrow(fromtlocs)){ text(fromtlocs[x,1],fromtlocs[x,2],pheno.col[nfrom(ctls)[x]]) }
  points(totlocs, cex=3, pch=20,col="white")
  points(totlocs, cex=3)
  for(x in 1:nrow(totlocs)){ text(totlocs[x,1],totlocs[x,2],nto(ctls)[x]) }
}

