#
# ctl.network.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Jan, 2012
# first written Oct, 2011
# 
# Network routines for CTL analysis
#

write.marker.attributes <- function(QTLscan, mapinfo){
  if(missing(QTLscan)) stop("argument 'CTLobject' is missing, with no default")
  if(!missing(mapinfo)){
    chr.edges    <- get.chr.edges(mapinfo)+.5
    chr          <- 1
    markernames  <- colnames(QTLscan)
  
    cat("", file="node_attributes.txt")
    #Connections between genetic markers
    for(m1 in 1:(length(markernames)-1)){
      cat(paste(markernames[m1],"\tCHR",chr,"\n",sep=""),file="node_attributes.txt",append=TRUE)
      if((((2*m1+1)/2) %in% chr.edges)){
        chr <- chr + 1
      }
    }
    #Annotate the last marker on the last chromosome
    cat(paste(markernames[length(markernames)],"\tCHR",chr,"\n",sep=""),file="node_attributes.txt",append=TRUE)
    cat("Wrote attributes for",length(markernames),"markers on",chr,"chromosomes\n")
  }else{
    cat("[WARNING] No nodeattributes.txt produced, Please supply a mapinfo file\n")
  }
}

QTLnetwork <- function(CTLobject, mapinfo, lod.threshold = 5, filename){
  if(missing(CTLobject)) stop("argument 'CTLobject' is missing, with no default")
  qtls <- ctl.qtlmatrix(CTLobject)
  if(!missing(mapinfo)){
    chr.edges   <- get.chr.edges(mapinfo)+.5
    write.marker.attributes(qtls, mapinfo)
  }
  if(missing(filename)){ 
    filename <- "network_qtl.sif" 
    cat("[WARNING] Using default filename:",filename,"\n")
  }
  idx         <- 1
  edges       <- 0
  ncons       <- 0
  chr         <- 1
  
  markernames <- colnames(qtls)
  traitnames  <- rownames(qtls)

  if(missing(filename)) 
  cat("",file=filename)
  #Print the CHR interconnections between genetic markers
  
  for(m1 in 1:(length(markernames)-1)){
    if(missing(mapinfo)){
      cat(paste(markernames[m1],"\tCHR\t",markernames[m1+1],"\t",1,"\n",sep=""),file=filename,append=TRUE)
      ncons <- ncons+1
    }else{
      if(!(((2*m1+1)/2) %in% chr.edges)){
        cat(paste(markernames[m1],"\tCHR",chr,"\t",markernames[m1+1],"\t",1,"\n",sep=""),file=filename,append=TRUE)
        ncons <- ncons+1
      }else{
        chr <- chr + 1
      }
    }
  }
  cat("Wrote",ncons," marker interconnections to",filename,"\n")
  #Print the QTLs between genetic markers and traits
  for(x in apply(qtls,1,function(x){which(x > lod.threshold)})){
    for(marker in names(x)){
      cat(paste(ctl.name(CTLobject[[idx]]),"\tQTL\t",marker,"\t",round(qtls[idx,marker],digits=1),"\n",sep=""),file=filename,append=TRUE)
      edges <- edges+1
    }
    idx <- idx + 1
  }
  cat("Wrote",edges,"QTL edges to",filename,"\n")
}

CTLnetwork <- function(CTLobject, mapinfo, lod.threshold = 5, add.qtl = TRUE, verbose = FALSE){
  if(missing(CTLobject)) stop("argument 'CTLobject' is missing, with no default")
  edges        <- 0
  edge_p_count <- NULL
  edge_m_count <- NULL

  cat("",file="network_full.sif")
  cat("",file="network_summary.sif")
  
  for(CTL in CTLobject){
    for(x in 1:ncol(CTL$ctl)){
      if(is.null(CTL$p)) stop("No permutations found, need permutations to determine likelihood\n")
      for(id in which(abs(CTL$ctl[,x]) > lod.threshold)){
        edge_name <- paste(ctl.name(CTL),"CTL",ctl.names(CTLobject)[x],sep="\t")
        if(edge_name %in% edge_p_count[,1]){
          edgeid <- which(edge_p_count[,1] %in% edge_name)
          if(verbose) cat("Duplicate edge",edgeid,":", edge_name,"\n")
          edge_p_count[edgeid,2] <- as.numeric(edge_p_count[edgeid,2])+abs(CTL$ctl[id,x])
          edge_p_count[edgeid,3] <- as.numeric(edge_p_count[edgeid,3])+1
        }else{
          edge_p_count <- rbind(edge_p_count,c(edge_name, CTL$ctl[id,x], 1))
        }
        cat(ctl.name(CTL),"CTL",rownames(CTL$ctl)[id],CTL$ctl[id,x],colnames(CTL$ctl)[x],"\n",file="network_full.sif",append=TRUE)
        edges <- edges+1
      }
    }
    if(verbose) cat("Analysed phenotype:",attr(CTL,"name"),"\n")
  }
  cat("Wrote",edges,"CTL edges to network_full.sif\n")
  if(!is.null(edge_p_count)){
    
    for(x in 1:nrow(edge_p_count)){
      cat(edge_p_count[x,1],"\t",round(as.numeric(edge_p_count[x,2]),digits=1),"\t",edge_p_count[x,3],"\n",file="network_summary.sif",append=TRUE,sep="")
    }
    cat("Wrote",nrow(edge_p_count),"CTL edges to network_summary.sif\n")
  }
  if(add.qtl) QTLnetwork(CTLobject, mapinfo, lod.threshold, "network_full.sif")
  if(add.qtl) QTLnetwork(CTLobject, mapinfo, lod.threshold, "network_summary.sif")
  invisible(unique(edge_p_count))
}

nodesToGenes <- function(nodenames, spotAnnotation){
  if(missing(nodenames)) stop("No names  to get annotation for supplied")
  if(missing(spotAnnotation)) stop("No annotation file supplied")
  cnt <- 0
  not_cnt <- 0
  genenames <- NULL
  for(x in nodenames){
    if(x %in% spotAnnotation$SPOT_ID){
      id <- which(spotAnnotation$SPOT_ID %in% x)
      genenames <- c(genenames,as.character(spotAnnotation$ORF[id[1]]))
      cnt <- cnt+1
    }else{
      not_cnt <- not_cnt+1
    }
  }
  cat("Found annotation for",cnt,"probes\n")
  cat("Not found",not_cnt,"probes\n")
  genenames
}

#Maps probeannot$SPOT_ID to the nodenames
write.node.attributes <- function(nodenames, spotAnnotation){
  if(missing(nodenames)) stop("No names  to get annotation for supplied")
  if(missing(spotAnnotation)) stop("No annotation file supplied")
  cnt <- 0
  cat("",file="node.attributes")
  for(x in nodenames){
    if(x %in% spotAnnotation$SPOT_ID)
    id <- which(spotAnnotation$SPOT_ID %in% x)
    desc_split <- strsplit(strsplit(as.character(spotAnnotation$DESCRIPTION[id]),"|",fixed=TRUE)[[1]],"=")
    if(length(desc_split) ==3){
      cat(x,"\t",as.character(spotAnnotation$ORF[id[1]]),"\t",desc_split[[1]][2],"\t",desc_split[[2]][2],"\t",desc_split[[3]][2],"\n",file="node.attributes",append=TRUE,sep="")
    }else{
      cat(x,"\t",as.character(spotAnnotation$ORF[id[1]]),"\tUnknown\tUnknown\tUnknown\n",file="node.attributes",append=TRUE,sep="")
    }
    cnt <- cnt+1
  }
  cat("Wrote annotation for",length(nodenames),"probes to node.attributes\n")
}

# end of ctl.network.R
