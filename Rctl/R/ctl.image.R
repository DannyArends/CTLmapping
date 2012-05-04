#
# ctl.image.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Jan, 2012
# first written Nov, 2010
# 
# Image plot routines for CTL analysis
#

image.CTLobject <- function(x, against = c("markers","phenotypes"), onlySignificant = FALSE, significance = 0.05, col=c("white",gray.colors(10)[10:1]), do.grid=TRUE, grid.col = "black", verbose = FALSE, add=FALSE, ...){
  colorrange <- col
  mymatrix <- NULL
  mynames <- NULL
  for(p in 1:length(x)){
    if(onlySignificant){
      maxes <- apply(abs(x[[p]]$ctl),1,max)
      if(max(maxes) > getPermuteThresholds(x[[p]], significance)[1]){
        mymatrix <- rbind(mymatrix,apply(x[[p]]$l,2,sum))
        mynames <- c(mynames,attr(x[[p]]$ctl,"name"))  
      }
    }else{
      mymatrix <- rbind(mymatrix,apply(x[[p]]$l,2,sum))
      mynames <- c(mynames,attr(x[[p]]$ctl,"name"))
    }
  }
  rownames(mymatrix) <- mynames
  mainlabel <- paste("CTL phenotypes vs",against[1],"at P-value <",significance)
  internal.image(mymatrix, colorrange, mainlabel,do.grid, grid.col)
}

qtlimage <- function(CTLscan, do.grid = TRUE, grid.col = "black", verbose = FALSE){
  colorrange <- c("white",gray.colors(10)[10:1])
  mymatrix <- NULL
  mynames <- NULL
  for(num in 1:length(CTLscan)){ 
    mymatrix <- rbind(mymatrix,as.numeric(unlist(CTLscan[[num]]$qtl))) 
    mynames <- c(mynames,attr(CTLscan[[num]]$ctl,"name"))  
  }
  rownames(mymatrix) <- mynames
  colnames(mymatrix) <- names(CTLscan[[1]]$qtl)
  internal.image(mymatrix, colorrange, "QTLs", do.grid, grid.col)
  invisible(mymatrix)
}

internal.image <- function(mymatrix, colorrange, mainlabel, do.grid, grid.col, add=FALSE){
  if(!is.null(mymatrix)){ 
    image(1:ncol(mymatrix),1:nrow(mymatrix),t(mymatrix), main=mainlabel, yaxt="n", 
          xaxt="n", ylab="", xlab="",col=c("white",gray.colors(4)[4:1]), cex.main=0.7, 
          breaks = c(0,2,4,8,10,100),add=add)
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
  mymatrix <- attr(x,"qtl")
  if(onlySignificant){
    mainlabel <- paste("QTL heatmap at P-value <", significance)
  }else{
    mainlabel <- paste("QTL heatmap")
  }
  internal.image(mymatrix, colorrange, mainlabel,do.grid, grid.col)
}

# end of ctl.image.R
