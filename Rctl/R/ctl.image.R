#
# ctl.image.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Nov, 2010
# 
# Image plot routines for CTL analysis
#

image.CTLobject <- function(x, against = c("markers","phenotypes"), significance = 0.05, col=whiteblack, do.grid=TRUE, grid.col = "white", verbose = FALSE, add=FALSE, breaks = c(0, 1, 2, 3, 10, 10000), marker_info, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  colorrange <- col
  mymatrix <- CTLprofiles(x, against, significance)
  mainlabel <- paste("CTL phenotypes vs",against[1],"at P-value <",significance)
  internal.image(mymatrix, colorrange, mainlabel,do.grid, grid.col, breaks=breaks, marker_info=marker_info)
}

qtlimage <- function(x, do.grid = TRUE, grid.col = "white", verbose = FALSE, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  colorrange <- c("white",gray.colors(10)[10:1])
  mymatrix <- NULL
  mynames <- NULL
  for(num in 1:length(x)){ 
    mymatrix <- rbind(mymatrix,as.numeric(unlist(x[[num]]$qtl))) 
    mynames <- c(mynames,ctl.name(x[[num]]))  
  }
  rownames(mymatrix) <- mynames
  colnames(mymatrix) <- names(x[[1]]$qtl)
  internal.image(mymatrix, colorrange, "QTLs", do.grid, grid.col, ...)
  invisible(mymatrix)
}

internal.image <- function(mymatrix, colorrange = whiteblack, mainlabel="Image", do.grid = TRUE, grid.col = 'white', add=FALSE, breaks = c(0, 1, 2, 3, 10, 10000), marker_info){
  if(missing(mymatrix)) stop("argument 'mymatrix' is missing, with no default")
  if(!is.null(mymatrix)){ 
    image(1:ncol(mymatrix),1:nrow(mymatrix),t(mymatrix), main=mainlabel, yaxt="n", 
          xaxt="n", ylab="", xlab="",col=c("white",gray.colors(4)[4:1]), cex.main=0.7, 
          breaks = breaks,add=add)
    axis(2,rownames(mymatrix),at=1:nrow(mymatrix),las=2,cex.axis=0.5)
    if(do.grid){
      abline(h=seq(-0.5,nrow(mymatrix)+0.5,1),col=grid.col,lwd=1)
      if(missing(marker_info)) abline(v=seq(-0.5,ncol(mymatrix)+0.5,1),col=grid.col,lwd=1)
    }
    if(missing(marker_info)){
      axis(1,colnames(mymatrix),at=1:ncol(mymatrix),las=2,cex.axis=0.6)
    }else{
      addChromosomeLines(marker_info)
    }
    box()
  }
  invisible(mymatrix)
}

addChromosomeLines <- function(markerinfo, col='black'){
  chr_start <- c()
  chr_ends <- c()
  for(x in unique(markerinfo[,1])){
    chr_start <- c(chr_start,min(which(markerinfo[,1]==x)))
    chr_ends  <- c(chr_ends, max(which(markerinfo[,1]==x)))
  }
  abline(v=.5 + chr_ends, lwd=1, lty=2, col=col)
  for(x in 1:length(chr_start)){
    axis(1,at=(chr_start[x]+chr_ends[x])/2,paste("Chr",unique(markerinfo[,1])[x]),cex.axis=0.3)
  }
  box()
}

# end of ctl.image.R
