#
# ctl.plot.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Nov, 2010
# 
# Plotting routines for CTL analysis
#

CTLasLOD <- function(CTLscan, QTLscores, main, do.legend=TRUE){
  if(missing(CTLscan)) stop("argument 'CTLscan' is missing, with no default")
  if(missing(main)) main <- paste("Comparison CTL:QTL of",ctl.name(CTLscan))
  
  op <- par(mfrow=c(2,1))
  plot.CTLscan(CTLscan,do.legend=do.legend)
  CTLscores <- CTLtoLODvector(CTLscan)
  if(!missing(QTLscores)){
    QTLscores <- as.numeric(unlist(QTLscores))
    plot(c(0,length(CTLscores)),c(0,max(c(CTLscores,QTLscores))),type='n', main=main, ylab="LOD",xlab="Marker")
    points(CTLscores,type='l',col="black",lwd=3)
    points(QTLscores,type='l',col="red",lwd=2,lty=1)
    if(do.legend) legend("topleft",c("CTL","QTL"),col=c("black","red"),lty=c(1,1),lwd=c(3,2))
  }else{
    plot(c(0,length(CTLscores)),c(0,max(c(CTLscores))),type='n', main=main, ylab="LOD",xlab="Marker")
    points(CTLscores,type='l',col="black",lwd=3)
  }
}

plot.CTLobject <- function(x, pheno.col=1:length(x), ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  if(length(x) == 1){
    return(plot.CTLscan(x[[1]],...))
  }
  if(length(pheno.col) == 1){
    return(plot.CTLscan(x[[pheno.col]],...))
  }
  return(image.CTLobject(x,...))
}

plotExpression <- function(genotypes, phenotypes, traits=c("X3.Hydroxypropyl", "X3.Methylthiopropyl"), markers=1, method = c("pearson", "kendall", "spearman"), do.plot = TRUE, verbose = FALSE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  if(length(traits) != 2) stop("argument 'traits' needs to be of length 2")
  
  ids <- c(which(colnames(phenotypes)==traits[1]),which(colnames(phenotypes)==traits[2]))
  result <- NULL
  
  for(m in markers){
    if(do.plot && length(markers)==1) plot(x=phenotypes[,ids[1]],y=phenotypes[,ids[2]],col=genotypes[,m],pch=20,xlab=traits[1],ylab=traits[2])
    idx1 <- which(genotypes[,m]==1)
    idx2 <- which(genotypes[,m]==2)
    a   <- cor(phenotypes[idx1,ids[1]],phenotypes[idx1,ids[2]],use="pair",method = method[1])
    b   <- cor(phenotypes[idx2,ids[2]],phenotypes[idx2,ids[1]],use="pair",method = method[1])
    ctl <- (a - b)^2
    if(verbose) cat("Correlation: ",a,"\t",b,"\tCTL: ",ctl,"\n")
    result <- rbind(result, c(a,b,ctl))
  }
  
  if(length(markers)>1 && do.plot){
    plot(x=c(0,length(markers)),y=c(-1,1),main=paste("CTL ",traits[2]),ylab="Correlation",xlab="Marker",t='n')
    points(result[,1],col="green",t='l')
    points(result[,2],col="red",t='l')
    points(result[,3],col="black",t='l',lwd=2)
  }
  return(result)
}

plot.CTLscan2 <- function(x, addQTL = TRUE, onlySignificant = TRUE, significance = 0.05, do.legend=TRUE, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  if(!is.null(x$perms) && !is.nan(getPermuteThresholds(x,significance)[1])){
    mysign <- as.numeric(which(apply(abs(x$ctl),1,max) > getPermuteThresholds(x,significance)[1]))
  }else{
    mysign <- as.numeric(which(apply(abs(x$ctl),1,max) > -log10(significance)))
  }
  if(length(mysign) ==0 || onlySignificant == FALSE){
    mysign <- 1:nrow(x$ctl)
    do.legend=FALSE
  }
  if(is.null(x$ctl)){
    CTLmatrix <- matrix(abs(x$ctl[mysign, ]),length(mysign),ncol(x$ctl))
  }else{
    CTLmatrix <- matrix(sign(x$ctl[mysign, ]) * x$ctl[mysign, ],length(mysign),ncol(x$ctl))
  }
  summarized <- apply(CTLmatrix,2,sum)
  plot(c(0,ncol(x$ctl)),c(min(c(summarized,x$qtl)),max(c(summarized,x$qtl))), type='n',xlab="Marker", ylab="-log10(P-value)", main=paste("Phenotype contribution to CTL of",ctl.name(x)),...)
  p <- rep(0,ncol(x$ctl))
  i <- 1;
  mycolors <- rep("black",nrow(CTLmatrix))
  apply(CTLmatrix,1,
    function(d){
      points(d,type='b',col=mycolors[i],pch=i,lwd=2)
      p <<- p + d
      i <<- i + 1
    }
  )
  if(do.legend){
    legend("topleft",c("QTL",paste("CTL",rownames(x$ctl)[mysign])),col=c("red",mycolors),lwd=2,pch=c(NA,1:length(mycolors)),cex=1.2)
  }
  if(!is.null(x$ctl)){
    n <- dim(x$ctl)[2]
    abline(h=-log10(c(0.05/n,0.01/n,0.001/n)),col=c("red","orange","green"),lty=2)
    abline(h=log10(c(0.05/n,0.01/n,0.001/n)),col=c("red","orange","green"),lty=2)
    legend("topright",as.character(paste("CTL-FDR:",c(0.05,0.01,0.001),"%")),col=c("red","orange","green"),lty=rep(2,3),lwd=1,cex=1.2)
  }
  points(as.numeric(x$qtl),type='l',lwd=2,col="red")
  rownames(CTLmatrix) <- rownames(x$ctl)[mysign]
  invisible(CTLmatrix)
}

plot.CTLscan <- function(x, addQTL = TRUE, onlySignificant = TRUE, significance = 0.05, do.legend=TRUE, cex.legend=1.0, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  mysign <- as.numeric(which(apply(abs(x$ctl),2,max) > -log10(significance)))
  if(length(mysign) ==0 || onlySignificant == FALSE){
    mysign <- 1:ncol(x$ctl)
    do.legend=FALSE
  }
  CTLmatrix <- matrix(x$ctl[,mysign],nrow(x$ctl),length(mysign))
  summarized <- apply(CTLmatrix,1,sum)
  x$qtl[is.infinite(x$qtl)] <- max(x$qtl[is.finite(x$qtl)])
  cat("MAX:",max(c(5,summarized)),"\n")
  plot(c(0,nrow(x$ctl)),c(min(c(summarized,x$qtl)),max(c(5,summarized,x$qtl))), type='n',xlab="Marker", ylab="-log10(P-value)", main=paste("Phenotype contribution to CTL of",ctl.name(x)), ...)
  p <- rep(0,nrow(x$ctl))
  i <- 1;
  mycolors <- topo.colors(ncol(CTLmatrix))
  apply(CTLmatrix,2,
    function(d){
     for(idx in 1:length(d)){
        rect(idx-0.5,p[idx],idx+0.5,p[idx]+d[idx],col=mycolors[i],lwd=0,lty=0)
      }
      p <<- p + d
      i <<- i + 1
    }
  )
  if(!is.null(x$ctl)){
    abline(h=-log10(c(0.05,0.01,0.001)),col=c("red","orange","green"),lty=2)
    legend("topright",as.character(paste("CTL-FDR:",c(0.05,0.01,0.001),"%")),col=c("red","orange","green"),lty=rep(2,3),lwd=1,cex=cex.legend)
  }
  if(do.legend){
    legend("topleft",colnames(x$ctl)[mysign],col=mycolors,lwd=1,cex=cex.legend)
  }
  points(summarized,type='l',lwd=1)
  points(as.numeric(x$qtl),type='l',lwd=2,col="red")
  colnames(CTLmatrix) <- colnames(x$ctl)[mysign]
  invisible(CTLmatrix)
}

chr_length <- function(map_info, chr = 1){ max(map_info[which(map_info[,1]==chr),2]) }

chr_total_length <- function(map_info, gap = 25){
  l <- 0
  for(x in unique(map_info[,1])){ l <- l + chr_length(map_info,x) + gap }
  (l-gap) #Gaps are between, so we don't need the last gap
}

a_loc <- function(map_info, id=1, gap = 25){
  res <- NULL
  for(x in 1:nrow(map_info)){ res <- c(res, m_loc(map_info,x)) }
  res
}

m_loc <- function(map_info, id=1, gap = 25){
  chr <- map_info[id,1]
  l <- 0
  while((chr-1) > 0){
   l <- l + chr_length(map_info, (chr-1))+gap
   chr <- chr - 1
  }
  l + map_info[id,2]
}

plot.CTLscan3 <- function(x, map_info, ...){
  summarized <- apply(x$ctl,1,sum)
  plot(c(200, chr_total_length(map_info)),c(0, max(x$ctl,x$qtl)), type='n',xlab="",ylab="")
  loc <- NULL
  for(y in 1:nrow(map_info)){loc <- c(loc,m_loc(map_info,y))}

  for(y in c(3,4,5)){
    onchr <- which(map_info[,1]==y)
    points(loc[onchr],as.numeric(x$qtl)[onchr],type='l', lwd=3, lty=1)
    for(g in 1:ncol(x$ctl)){
      points(loc[onchr],as.numeric(x$ctl[onchr,g]),type='l',lwd=3, lty=4, col=(g+1))
    }
  }
  traitnamez <- gsub(".Mean","", colnames(x$ctl))
  legend("topleft",c(paste("QTL",gsub(".Mean","",ctl.name(x))),paste("CTL",traitnamez[1]),paste("CTL",traitnamez[2]),paste("CTL",traitnamez[3])),lwd=c(2, 2, 2, 2), lty=c(1,4,4,4),col=c("black",2,3,4))
}

plot.CTLpermute <- function(x, type="s", ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  #plot(seq(0,0.9,0.01),CTLscoretoPvalue(seq(0,0.9,0.01),x),main="CTL to P.value",xlab="CTL",ylab="Pvalue", type=type,...)
  significant <- print.CTLpermute(x)
  mycolors <- c("red","orange","green")
  idx <- 1
  for(y in c(.05,.01,.001)){
    lines(rbind(c(-1,y),c(significant[idx],y)),lty=2,col=mycolors[idx])
    lines(rbind(c(significant[idx],-1),c(significant[idx],y)),lty=2,col=mycolors[idx])
    idx <- idx+1
  }
  legend("topright",c("CTL-FDR: 5%","CTL-FDR: 1%","CTL-FDR: 0.1%"), col=mycolors, lty=2, lwd=1, cex=0.7)
}

# end of ctl.plot.R
