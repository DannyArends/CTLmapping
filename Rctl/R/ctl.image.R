#
# ctl.image.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Nov, 2010
# 
# Image plot routines for CTL analysis
#

image.CTLobject <- function(x, marker_info, against = c("markers","phenotypes"), significance = 0.05, col=whiteblack, do.grid=TRUE, grid.col = "white", verbose = FALSE, add=FALSE, breaks = c(0, 1, 2, 3, 10, 10000), ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  colorrange <- col
  mymatrix <- CTLprofiles(x, against, significance)
  main <- paste("CTL phenotypes vs",against[1],"at P-value <",significance)
  internal.image(mymatrix, marker_info = marker_info, colorrange, main = main, do.grid, grid.col, breaks=breaks)
}

qtlimage <- function(x, marker_info, do.grid = TRUE, grid.col = "white", verbose = FALSE, ...){
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
  internal.image(mymatrix, marker_info = marker_info, colorrange, "QTLs", do.grid = do.grid, grid.col = grid.col, ...)
  invisible(mymatrix)
}

internal.image <- function(mymatrix, marker_info, colorrange = whiteblack, main = "Image", do.grid = TRUE, grid.col = 'white', add=FALSE, breaks = c(0, 1, 2, 3, 10, 10000), cex.axis = 1.0, cex.main = 1.0, las = 2){
  if(missing(mymatrix)) stop("argument 'mymatrix' is missing, with no default")
  if(!is.null(mymatrix)){ 
    image(1:ncol(mymatrix),1:nrow(mymatrix),t(mymatrix), main = main, yaxt="n", 
          xaxt="n", ylab="", xlab="", col = c("white", gray.colors(4)[4:1]), cex.main=cex.main, 
          breaks = breaks, add = add)
    axis(2,rownames(mymatrix),at=1:nrow(mymatrix), las = las, cex.axis = cex.axis)
    if(do.grid){
      abline(h=seq(-0.5,nrow(mymatrix)+0.5,1),col=grid.col,lwd = 1)
      if(missing(marker_info)) abline(v=seq(-0.5,ncol(mymatrix)+0.5,1),col=grid.col,lwd=1)
    }
    if(missing(marker_info)){
      axis(1, colnames(mymatrix),at=1:ncol(mymatrix), las = las, cex.axis = cex.axis)
    }else{
      addChromosomeLines(marker_info)
    }
    box()
  }
  invisible(mymatrix)
}

addChromosomeLines <- function(markerinfo, col='black', cex.axis = 1.0){
  chr_start <- c()
  chr_ends <- c()
  for(x in unique(markerinfo[,1])) {
    chr_start <- c(chr_start,min(which(markerinfo[,1]==x)))
    chr_ends  <- c(chr_ends, max(which(markerinfo[,1]==x)))
  }
  abline(v=.5 + chr_ends, lwd=1, lty=2, col=col)
  for(x in 1:length(chr_start)){
    axis(1, at=(chr_start[x]+chr_ends[x])/2, paste("Chr",unique(markerinfo[,1])[x]), cex.axis = cex.axis)
  }
  box()
}

# end of ctl.image.R
