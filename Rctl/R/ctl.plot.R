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

plotExpression <- function(genotypes, phenotypes, traits=c("X3.Hydroxypropyl", "X3.Methylthiopropyl"),
          markers=1, method = c("pearson", "kendall", "spearman"), do.plot = TRUE, verbose = FALSE){

  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  if(length(traits) != 2) stop("argument 'traits' needs to be of length 2")
  
  ids <- c(which(colnames(phenotypes)==traits[1]),which(colnames(phenotypes)==traits[2]))
  result <- NULL
  
  for(m in markers){
    if(do.plot && length(markers)==1){
      plot(x=phenotypes[,ids[1]],y=phenotypes[,ids[2]],col=genotypes[,m],pch=20,xlab=traits[1],ylab=traits[2])
    }
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

plot.CTLscan <- function(x, mapinfo = NULL, type = c("barplot","gwas","line"), onlySignificant = TRUE, significance = 0.05, gap = 25, plot.cutoff = FALSE, do.legend=TRUE, cex.legend=1.0, ydim=NULL, ylab="-log10(P-value)", ...){

  if(missing(x) || is.null(x)) stop("argument 'x' is missing, with no default")

  significant <- as.numeric(which(apply(abs(x$ctl), 2, max) > -log10(significance)))

  if(length(significant) ==0 || onlySignificant == FALSE){ significant <- 1:ncol(x$ctl); do.legend = FALSE }

  if(is.null(mapinfo)){
    maxX <- nrow(x$ctl) ; pointsx <- 1:nrow(x$ctl)
  }else{ # We have a mapinfo
    maxX <- chr_total_length(mapinfo) ; pointsx <- a_loc(mapinfo) ;
    pointsx <- c(pointsx, max(pointsx)) # Note: Add the max to pointsx
  }
  
  ctlsubset <- matrix(x$ctl[, significant], nrow(x$ctl), length(significant))
  colnames(ctlsubset) <- colnames(x$ctl)[significant]

  summarized <- apply(ctlsubset, 1, sum) # Summarized scores for all significant
  x$qtl[is.infinite(x$qtl)] <- max(x$qtl[is.finite(x$qtl)])
  if(is.null(ydim)){ 
    maxy <- max(c(7.5, ctlsubset, x$qtl))
    if(type[1] == "barplot") maxy <- max(c(7.5, summarized, x$qtl)) # Maximum of barplot is summarized
    ydim <- c(-maxy, maxy)
  }
  plot(c(0.5, maxX+0.5), ydim, type='n',xlab="", ylab=ylab, ...)
  points(pointsx, rep(0, length(pointsx)), lwd = 1, pch="|",cex = 0.2)

  i <- 1;
  ntraits  <- ncol(ctlsubset)
  nmarkers <- nrow(ctlsubset)
  ltype    <- 'l'
  if(type[1]=="gwas") ltype <- 'h' # GWAS plot uses bars
  colfunc <- colorRampPalette(c("red", "blue", "darkgreen", "orange"))
  mycolors <- colfunc(ncol(ctlsubset))

  p  <- rep(0,nrow(ctlsubset))
  mx <- 0
  apply(ctlsubset, 2, function(d){
     if(type[1]=="barplot"){        # Summarized bar plot
       for(idx in 1:length(d)){
          if(is.null(mapinfo)){
            mx   <- idx; lbar <- 0.5; rbar <- 0.5;
          }else{
            mp <- mx; mx <- pointsx[idx]
            lbar <- abs(mp - mx)/2; rbar <- abs(pointsx[idx+1] - mx)/2;
          }
          rect(mx-lbar, -p[idx],mx+rbar, -(p[idx]+d[idx]), col=mycolors[i], lwd = 0, lty = 0)
        }
      }else{ # Line or GWAS plot
        if(is.null(mapinfo)){ # Without mapinfo coordinated match
          points(pointsx, -d, type = ltype, lwd = 1, col = mycolors[i])
        }else{ # Go through the chromosomes
          for(chr in unique(mapinfo[, 1])){
            idxes <- which(mapinfo[, 1] == chr)
            if(type[1]=="gwas"){  # Plot Gwas, chromosome colors alternate
              lty   <- 1
              col <- c("black","orange")[(as.numeric(chr) %% 2)+1]
            }else{
              lty <- i ; col <- mycolors[i]
            }
            points(pointsx[idxes], -d[idxes], type = ltype, lwd = 2, lty = lty, col = col)
          }
        }
      }
      p <<- p + d
      i <<- i + 1
    } )

  # Plot the cut-off line at -log10(significance)
  if(plot.cutoff){
    abline(h=-log10(c(significance / nmarkers)), col=c("green"), lty=2)
    abline(h= log10(c(significance)), col=c("green"), lty=2)
    mleg <- as.character(paste("FDR:",c(significance),"%"))
  }
  # Plot the legend(s)
  if(do.legend){
    if(plot.cutoff) legend("topright", mleg, col=c("green"), lty=rep(2), lwd=1, cex=cex.legend, bty='n')
    lty <- 1:ntraits
    if(type[1] == "barplot") lty <- 1
    legend("topleft", colnames(ctlsubset), col=mycolors, lwd=1, lty=lty, cex=cex.legend, bty='n')
  }
  # Plot the summarized lines and QTLs
  if(is.null(mapinfo)){
    if(type[1] == "barplot") points(pointsx, -summarized,type='l',lwd=1)
    points(pointsx, as.numeric(x$qtl), type=ltype,lwd=2, col="black")
  }else{ # With mapinfo, plot chromosomes and markers
    for(chr in unique(mapinfo[, 1])){
      idxes <- which(mapinfo[, 1] == chr)
      if(type[1] == "barplot") points(pointsx[idxes], -summarized[idxes],type = ltype,lwd=1)
      points(pointsx[idxes], as.numeric(x$qtl)[idxes], type = ltype, lwd=2, col="black")
    }
  }
  invisible(ctlsubset)
}

chr_length <- function(mapinfo, chr = 1){ 
  max(mapinfo[which(mapinfo[, 1]==chr), 2]) 
}

chr_total_length <- function(mapinfo, gap = 25){
  l <- 0
  for(x in unique(mapinfo[, 1])){ l <- l + chr_length(mapinfo,x) + gap }
  (l-gap) #Gaps are between, so we don't need the last gap
}

a_loc <- function(mapinfo, gap = 25){
  res <- NULL
  for(x in 1:nrow(mapinfo)){ res <- c(res, m_loc(mapinfo,x, gap = gap)) }
  res
}

m_loc <- function(mapinfo, id = 1, gap = 25){
  chr <- as.numeric(mapinfo[id, 1])
  l <- 0
  while((chr-1) > 0){
   l <- l + chr_length(mapinfo, (chr-1)) + gap
   chr <- chr - 1
  }
  l + mapinfo[id, 2]
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
  legend("topleft", c(paste("QTL", gsub(".Mean","",ctl.name(x))),paste("CTL",traitnamez[1]),paste("CTL",traitnamez[2]),paste("CTL",traitnamez[3])),lwd=c(2, 2, 2, 2), lty=c(1,4,4,4),col=c("black",2,3,4))
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
