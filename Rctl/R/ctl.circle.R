#
# ctl.circle.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends and Ritsert C. Jansen
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

draw.element <- function(x, y, title, cex=1, bg.col = "white", border.col="black"){
  points(cbind(x, y), cex=cex*4, pch=19, col=bg.col)
  points(cbind(x, y), cex=cex*4, col= border.col)
  text(x, y, title, cex=cex)
}

circle.loc <- function(nt, size = 1.0){
  locs <- matrix(nrow = nt, ncol = 2)
  phi  <- seq(0, 2 * pi, length = (nt+1))
  complex.circle <- complex(modulus = 1, argument = phi)
  for(j in 1:nt){ locs[j, ] <- c(Im(complex.circle[j]), Re(complex.circle[j])); }
  invisible(locs * size)
}

nfrom <- function(ctls){ return(unique(ctls[,1])) }
nto   <- function(ctls){ return(unique(ctls[,3])) }

mapinfotomarkerlocs <- function(mapinfo, gap, type=c("line","circle")){
  if(missing(mapinfo) || is.null(mapinfo)) stop("argument 'mapinfo' is missing, with no default")
  if(length(type) > 1) type = type[1]
  n.chr       <- unique(mapinfo[,1])
  chr.lengths <- NULL
  markerlocs  <- NULL
  for(chr in 1:length(n.chr)){ #Absolute length of the chromosomes
    ll <- mapinfo[lapply(unique(mapinfo[,1]),function(x){which(x==mapinfo[,1])})[[chr]],]
    chr.lengths <- c(chr.lengths, max(as.numeric(ll[,2]))-min(as.numeric(ll[,2])))
  }
  total.l <- ceiling(sum(chr.lengths) + (gap*length(n.chr)))
  if(type == "line"){
    cmmap   <- cbind(seq(1,total.l),rep(0,total.l))
  }else if(type=="circle"){
    cmmap   <- circle.loc(total.l+5, 0.7) # TODO: Bug - Markers at -cM chromosome 1 -> negative idx
  }else{ stop("Type not supported, (Options: line & circle)"); }  
  for(x in 1:nrow(mapinfo)){
    m.chr <- which(unique(mapinfo[,1]) %in% mapinfo[x,1])
    m.loc <- (ceiling(mapinfo[x,2])+1)
    if(m.chr > 1){ m.loc <- m.loc + sum(chr.lengths[1:(m.chr-1)]) + ((m.chr-1)*gap); }
    markerlocs <- rbind(markerlocs, cmmap[ceiling(m.loc),])
  }
  invisible(markerlocs)
}

ctl.circle <- function(CTLobject, mapinfo, pheno.col, significance=0.05, gap=50, cex=1, verbose=FALSE){
  if(missing(CTLobject) || is.null(CTLobject)) stop("argument 'CTLobject' is missing, with no default")
  if(missing(pheno.col)) pheno.col <- 1:length(CTLobject) 

  CTLobject  <- CTLobject[pheno.col]  #Scale down to pheno.col as input
  n.markers  <- nrow(CTLobject[[1]]$ctl)
  ctls       <- CTLnetwork(CTLobject, mapinfo, significance, verbose = verbose)

  if(is.null(ctls)) return() # No ctls found
  markerlocs <- circle.loc(n.markers, 0.7)
  if(!missing(mapinfo)) markerlocs <- mapinfotomarkerlocs(mapinfo, gap, "circle")
  fromtlocs  <- circle.loc(length(nfrom(ctls)), 1.0)
  totlocs    <- circle.loc(length(nto(ctls)), 0.4)
  plot(c(-1.1, 1.1), c(-1.1, 1.1), type = "n", axes = FALSE, xlab = "", ylab = "")
  points(markerlocs, pch=20, cex=(cex/2))   # Plot the markers
  for(x in 1:nrow(ctls)){                   # Plot the ctls
    from      <- fromtlocs[which(nfrom(ctls) %in% ctls[x,1]),]
    to   <- totlocs[which(nto(ctls) %in% ctls[x,3]),]
    via  <- markerlocs[ctls[x,2],]
    draw.spline(from, to, via, lwd=(ctls[x,4]/5)+1, col=ctls[x,1])
  } # All done now plot the trait elements
  for(x in 1:nrow(fromtlocs)){ 
    draw.element(fromtlocs[x,1], fromtlocs[x,2], pheno.col[nfrom(ctls)[x]], cex=cex) 
  }
  for(x in 1:nrow(totlocs)){ 
    draw.element(totlocs[x,1], totlocs[x,2], nto(ctls)[x], cex=cex) 
  }
}

