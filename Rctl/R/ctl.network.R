#
# ctl.network.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Dec, 2012
# first written Oct, 2011
# 
# Network routines for CTL analysis
#

plot.CTLnetwork <- function(x, main="Causal significance", ...){
  clrs <- rgb((1-x[,4]/max(x[,4])), (x[,4]+1)/(max(x[,4])+1), 0)
  plot(x[,6], x[,7], col=clrs, pch=20, xlab="LOD Trait 1", ylab="LOD Trait 2", main = main, ...)
}

CTLnetwork <- function(CTLobject, mapinfo, significance = 0.05, LODdrop = 2, what = c("names","ids"), short = FALSE, add.qtls = FALSE, file = "", verbose = TRUE){
  if(missing(CTLobject) || is.null(CTLobject)) stop("argument 'CTLobject' is missing, with no default")
  if(any(class(CTLobject)=="CTLscan")) CTLobject = list(CTLobject)
  if(length(what) > 1) what = what[1]

  results <- NULL
  significant <- CTLsignificant(CTLobject, significance, what = "ids")
  if(!is.null(significant)){
    all_m <- NULL; all_p <- NULL;
    nodefile=""; netfile = "";
    if(file != ""){
      netfile <- paste("ctlnet",file,".sif",sep="")
      nodefile <- paste("ctlnet",file,".nodes",sep="")
    }
    cat("",file=netfile); cat("",file=nodefile);
    if(verbose) cat("NETWORK.SIF\n")
    edges <- NULL
    for(x in 1:nrow(significant)){
      data    <- as.numeric(significant[x,])
      CTLscan <- CTLobject[[data[1]]]
      markern <- rownames(CTLscan$dcor)
      traitsn <- colnames(CTLscan$dcor)
      name    <- ctl.name(CTLscan)
      if(what=="ids"){
        tid     <- which(traitsn %in% ctl.name(CTLobject[[data[1]]]))
        name    <- paste("P",tid,sep="")
        markern <- paste("M",1:nrow(CTLobject[[data[1]]]$dcor), sep="")
        traitsn <- paste("P", 1:ncol(CTLobject[[data[1]]]$dcor), sep="")
      }
      if(add.qtls){ # Add QTL to the output SIF
        bfc    <- length(CTLscan$qtl)
        above  <- which(CTLscan$qtl > -log10(significance))
        qtlnms <- names(above); qtlmid <- 1
        for(m in above){
          cat(name,"\tQTL\t",markern[m],"\tQTL\t",CTLscan$qtl[m],"\n",sep="",file=netfile,append=TRUE)
          all_m  <- CTLnetwork.addmarker(all_m, mapinfo, markern[data[2]], qtlnms[qtlmid])
          qtlmid <- qtlmid+1
        }
      }
      lod <- CTLscan$ctl[data[2],data[3]]
      qlod1    <- CTLscan$qtl[data[2]]
      qlod2    <- qlod1
      edgetype <- NA
      if(length(CTLobject) >= data[3]){  # Edge type based on QTL LOD scores
        qlod2 <-CTLobject[[data[3]]]$qtl[data[2]]
        if((qlod1-qlod2) > LODdrop){
          edgetype <- 1
        }else if((qlod1-qlod2) < -LODdrop){
          edgetype <- -1
        }else{ edgetype <- 0; }
      }else{ cat("Warning: Phenotype", data[3], "from", data[1], "no CTL/QTL information"); }
      #Store the results
      results <- rbind(results, c(data[1], data[2], data[3], lod, edgetype, qlod1, qlod2))

      if(nodefile == "" && !verbose){ }else{
        if(short){
          edge <- paste(name,traitsn[data[3]])
          edgeI <- paste(traitsn[data[3]],name)
          if(!edge %in% edges && !edgeI %in% edges){
            cat(name, "\t", markern[data[2]],"\t", traitsn[data[3]],"\n",file=netfile, append=TRUE,sep="")
            edges <- c(edges,edge)
          }
        }else{
          cat(name, "\t", "CTL_", data[1],"_",data[3], "\t", markern[data[2]],file=netfile, append=TRUE,sep="")
          cat("\tCTL\t", lod, "\n", file=netfile, append=TRUE,sep="")
          cat(markern[data[2]], "\t", "CTL_", data[1],"_",data[3], "\t",file=netfile, append=TRUE,sep="")
          cat(traitsn[data[3]],"\tCTL\t", lod, "\n", file=netfile,append=TRUE,sep="")
        }
      }
      all_m <- CTLnetwork.addmarker(all_m, mapinfo, markern[data[2]], rownames(CTLscan$dcor)[data[2]])
      all_p <- unique(c(all_p, name, traitsn[data[3]]))
    }
    colnames(results) <- c("TRAIT1","MARKER","TRAIT2","LOD_C","CAUSAL","LOD_T1","LOD_T2")
    if(verbose) cat("NODE.DESCRIPTION\n")
    if(nodefile == "" && !verbose){ }else{
      for(m in all_m){ cat(m,"\n",    sep="", file=nodefile, append=TRUE); }
      for(p in all_p){ cat(p,"\tPHENOTYPE\n", sep="", file=nodefile, append=TRUE); }
    }
  }
  class(results) <- c(class(results),"CTLnetwork")
  invisible(results)
}

CTLnetwork.addmarker <- function(markers, mapinfo, name, realname){
  if(!missing(mapinfo)){
    id      <- which(rownames(mapinfo) %in% realname)
    fname   <- paste(name,"\tMARKER\t",mapinfo[id,1],"\t",mapinfo[id,2],sep="")
    markers <- unique(c(markers, fname))
  }
  return(markers)
}

postscript("test.eps", width = 19.0, height = 19.0, paper="special", horizontal=FALSE)
op <- par(mfrow=c(3,3))
for(x in 1:8){  ctl.lineplot(ctls, map_info, x, 0.0001,cex=2, col="darkorange") }
plot(1:10,1:10,t='n',xaxt='n',xlab="",,yaxt='n',ylab="")
for(x in 1:9){
  draw.element(3,10-(x - .5),as.character(x),cex=2, bg.col="lightgray")
  text(6,10-(x - .5),ctl.names(ctls)[x],cex=2)
}
dev.off()
# end of ctl.network.R
