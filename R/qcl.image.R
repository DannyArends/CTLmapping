#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Nov, 2010
# 
# Image plot routines for QCL analysis
#

image.QCLobject <- function(x, against = c("markers","phenotypes"), onlySignificant = FALSE, significance = 0.05, do.grid=TRUE, grid.col = "black", verbose = FALSE, ...){
  colorrange <- c("white",gray.colors(10)[10:1])
  mymatrix <- NULL
  mynames <- NULL
  for(p in 1:length(x)){
    if(verbose) cat("Processing:",p,"from QCL to LOD\n")
    lod <- QCLtoLODvector(x[[p]], against)
    if(onlySignificant){
      if(max(lod) > -log10(significance)){
        mymatrix <- rbind(mymatrix,lod)
        mynames <- c(mynames,attr(x[[p]]$qcl,"name"))  
      }
    }else{
      mymatrix <- rbind(mymatrix,lod)
      mynames <- c(mynames,attr(x[[p]]$qcl,"name"))
    }
  }
  rownames(mymatrix) <- mynames
  mainlabel <- paste("QCL phenotypes vs",against[1],"at P-value <",significance)
  internal.image(mymatrix, colorrange, mainlabel,do.grid, grid.col)
}

internal.image <- function(mymatrix, colorrange, mainlabel, do.grid, grid.col){
  if(!is.null(mymatrix)){ 
    image(1:ncol(mymatrix),1:nrow(mymatrix),t(mymatrix),
          main=mainlabel,
          yaxt="n",xaxt="n",ylab="", xlab="",col=colorrange,cex.main=0.7,breaks = c(0,1,2,3,4,5,6,7,8,9,10,100))
    axis(2,rownames(mymatrix),at=1:nrow(mymatrix),las=2,cex.axis=0.5)
    axis(1,colnames(mymatrix),at=1:ncol(mymatrix),las=2,cex.axis=0.5)
    if(do.grid){
      abline(h=seq(-0.5,nrow(mymatrix)+0.5,1),col=grid.col,lwd=1)
      abline(v=seq(-0.5,ncol(mymatrix)+0.5,1),col=grid.col,lwd=1)
    }
    box()
  }
  invisible(mymatrix)
}

QTLimage <- function(x, onlySignificant = FALSE, significance = 0.05, do.grid=TRUE, grid.col = "black", verbose = FALSE, ...){
  colorrange <- c("white",gray.colors(10)[10:1])
  mymatrix <- NULL
  mynames <- NULL
  for(p in 1:length(x)){
    if(verbose) cat("Processing QTL",p,"from QCLobject\n")
    lod <- x[[p]]$qtl
    if(onlySignificant){
      if(max(QCLtoLODvector(x[[p]], "markers")) > -log10(significance)){
        mymatrix <- rbind(mymatrix,lod)
        mynames <- c(mynames,attr(x[[p]]$qcl,"name"))  
      }
    }else{
      mymatrix <- rbind(mymatrix,lod)
      mynames <- c(mynames,attr(x[[p]]$qcl,"name"))
    }
  }
  rownames(mymatrix) <- mynames
  if(onlySignificant){
    mainlabel <- paste("QTL heatmap at P-value <", significance)
  }else{
    mainlabel <- paste("QTL heatmap")
  }
  internal.image(mymatrix, colorrange, mainlabel,do.grid, grid.col)
}