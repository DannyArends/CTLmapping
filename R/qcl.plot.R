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
    main <- paste("Comparison QCL:QTL of",attr(QCLscan$qcl,"name"))
  }
  op <- par(mfrow=c(2,1))
  plot.QCLscan(QCLscan,do.legend=do.legend)
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
    plot.QCLscan(x[[1]],...)
  }else{
    image.QCLobject(x,...)
  }
}

plot.QCLscan <- function(x, onlySignificant = TRUE, qcl.threshold =0.6, do.legend=TRUE, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  if(!is.null(x$p)){
    mysign <- as.numeric(which(apply(abs(x$qcl),1,max) > getPermuteThresholds(x)[1]))
  }else{
    mysign <- as.numeric(which(apply(abs(x$qcl),1,max) > qcl.threshold))
  }
  if(length(mysign) ==0){
    mysign <- 1:nrow(x$qcl)
    do.legend=FALSE
  }
  if(is.null(x$l)){
    QCLmatrix <- matrix(abs(x$qcl[mysign, ]),length(mysign),ncol(x$qcl))
  }else{
    QCLmatrix <- matrix(x$l[mysign, ],length(mysign),ncol(x$l))
  }
  summarized <- apply(QCLmatrix,2,sum)
  plot(c(0,ncol(x$qcl)),c(0,max(summarized)), type='n',xlab="Marker", ylab="-log10(P-value)", main=paste("Phenotype contribution to QCL of",attr(x$qcl,"name")),...)
  p <- rep(0,ncol(x$qcl))
  i <- 1;
  mycolors <- topo.colors(nrow(QCLmatrix))
  apply(QCLmatrix,1,
    function(d){
     for(idx in 1:length(d)){
        rect(idx-0.5,p[idx],idx+0.5,p[idx]+d[idx],col=mycolors[i],lwd=0,lty=0)
      }
      p <<- p + d
      i <<- i + 1
    }
  )
  if(!is.null(x$l)){
    abline(h=-log10(c(0.05,0.01,0.001)),col=c("red","orange","green"),lty=2)
    legend("topright",as.character(paste("QCL-FDR:",c(0.05,0.01,0.001),"%")),col=c("red","orange","green"),lty=rep(2,3),lwd=1,cex=0.7)
  }
  if(do.legend){
    legend("topleft",rownames(x$qcl)[mysign],col=mycolors,lwd=1,cex=0.7)
  }
  points(summarized,type='l',lwd=1)
  rownames(QCLmatrix) <- rownames(x$qcl)[mysign]
  invisible(QCLmatrix)
}
