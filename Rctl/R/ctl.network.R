#
# ctl.network.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Oct, 2011
# 
# Network routines for CTL analysis
#

write.marker.attributes <- function(mapinfo){
  if(!missing(mapinfo)){
    chr.edges    <- get.chr.edges(mapinfo)+.5
    chr          <- 1
    markernames  <- rownames(mapinfo)
  
    #Connections between genetic markers
    for(m1 in 1:(length(markernames)-1)){
      cat(paste(markernames[m1],"\tCHR",chr,"\n",sep=""))
      if((((2*m1+1)/2) %in% chr.edges)){
        chr <- chr + 1
      }
    }
    #Annotate the last marker on the last chromosome
    cat(paste(markernames[length(markernames)],"\tCHR",chr,"\n",sep=""))
  }
}

CTLnetwork <- function(CTLobject, significance = 0.05, what = c("names","ids"), add.qtls = FALSE, mapinfo, file = ""){
  if(length(what) > 1) what = what[1]
  significant <- CTLsignificant(CTLobject, significance, what = "ids")
  if(!is.null(significant)){
    all_m <- NULL; all_p <- NULL;
    nodefile="";   netfile = "";
    if(file != ""){
      netfile <- paste("ctlnet",file,".sif",sep="")
      nodefile <- paste("ctlnet",file,".nodes",sep="")
    }
    cat("",file=netfile); cat("",file=nodefile);
    cat("NETWORK.SIF\n")
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
      if(add.qtls){
        bfc <- length(CTLscan$qtl)
        above <- which(CTLscan$qtl > -log10((0.05/bfc)))
        qtlmarkernames <- names(above); qtlmid <- 1
        for(m in above){
          cat(name, "\tQTL\t", markern[m],"\tQTL\t", CTLscan$qtl[m], "\n", sep="", file=netfile, append=TRUE)
          all_m <- CTLnetwork.addmarker(all_m, mapinfo, markern[data[2]], qtlmarkernames[qtlmid])
          qtlmid <- qtlmid+1
        }
      }
      cat(name, "\t", "CTL_", data[1],"_",data[3], "\t", markern[data[2]],"\tCTL\t", CTLscan$ctl[data[2],data[3]], "\n",sep="",file=netfile,append=TRUE)
      cat(markern[data[2]], "\t", "CTL_", data[1],"_",data[3], "\t", traitsn[data[3]],"\tCTL\t", CTLscan$ctl[data[2],data[3]], "\n",sep="",file=netfile,append=TRUE)
      all_m <- CTLnetwork.addmarker(all_m, mapinfo, markern[data[2]], rownames(CTLscan$dcor)[data[2]])
      all_p <- unique(c(all_p, name, traitsn[data[3]]))
    }
    cat("NODE.DESCRIPTION\n")
    for(m in all_m){ cat(m,"\n",    sep="", file=nodefile, append=TRUE); }
    for(p in all_p){ cat(p,"\tPHENOTYPE\n", sep="", file=nodefile, append=TRUE); }
  }
}

CTLnetwork.addmarker <- function(markers, mapinfo, name, realname){
  if(!missing(mapinfo)){
    id      <- which(rownames(mapinfo) %in% realname)
    fname   <- paste(name,"\tMARKER\t",mapinfo[id,1],"\t",mapinfo[id,2],sep="")
    markers <- unique(c(markers, fname))
  }
  return(markers)
}

# end of ctl.network.R
