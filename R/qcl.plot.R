#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Oct, 2011
# first written Nov, 2010
# 
# Plotting routines for QCL analysis
#

QCLasLOD <- function(QCLscan, QTLscores, main, do.legend=TRUE){
  if(missing(QCLscan)) stop("argument 'QCLscan' is missing, with no default")

  if(missing(main)){
    main <- paste("Comparison QCL:QTL of",attr(QCLscan$s,"name"))
  }
  QCLscores <- QCLtoLODvector(QCLscan)
  if(!missing(QTLscores)){
    plot(c(0,length(QCLscores)),c(0,max(c(QCLscores,QTLscores))),type='n', main=main, ylab="LOD",xlab="Marker")
    points(QCLscores,type='l',col="black",lwd=3)
    points(QTLscores,type='l',col="red",lwd=2,lty=1)
    if(do.legend) legend("topleft",c("QCL","QTL"),col=c("black","red"),lty=c(1,1),lwd=c(3,2))
  }else{
    plot(c(0,length(QCLscores)),c(0,max(c(QCLscores))),type='n', main=main, ylab="LOD",xlab="Marker")
    points(QCLscores,type='l',col="black",lwd=3)
  }
}

plot.QCLobject <- function(x, ...){
  if(length(x) == 1){
    plot.QCLscan(x[[1]])
  }else{
    image.QCLobject(x)
  }
}

plot.QCLscan <- function(x, onlySignificant = TRUE, do.legend=TRUE, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  mysign <- as.numeric(which(apply(abs(x$s),1,max) > getPermuteThresholds(x$p)[1]))
  if(length(mysign) ==0){
    mysign <- 1:nrow(x$s)
  }
  QCLmatrix <- matrix(x$l[mysign, ],length(mysign),ncol(x$l))
  summarized <- apply(QCLmatrix,2,sum)
  plot(c(0,ncol(x$s)),c(0,max(summarized)), type='n',xlab="Marker", ylab="-log10(P-value)", main=paste("Phenotype contribution to QCL of",attr(x$s,"name")),...)
  p <- rep(0,ncol(x$l))
  i <- 1;
  mycolors <- terrain.colors(nrow(QCLmatrix))
  apply(QCLmatrix,1,
    function(d){
     for(idx in 1:length(d)){
        rect(idx-0.5,p[idx],idx+0.5,p[idx]+d[idx],col=mycolors[i])
      }
      p <<- p + d
      i <<- i + 1
    }
  )
  if(do.legend) legend("topleft",rownames(x$l)[mysign],col=mycolors,lwd=1,cex=0.7)
  points(summarized,type='l',lwd=2)
  invisible(summarized)
}
