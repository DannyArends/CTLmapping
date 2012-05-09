#
# ctl.network.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Jan, 2012
# first written Oct, 2011
# 
# Network routines for CTL analysis
#

write.marker.attributes <- function(QTLscan, chr.edges){
  chr   <- 1
  namez <- colnames(QTLscan)
  cat("", file="node_attributes.txt")
  #Connections between genetic markers
  for(m1 in 1:(length(namez)-1)){
    if(!missing(chr.edges)){
      cat(paste(namez[m1],"\tCHR",chr,"\n",sep=""),file="node_attributes.txt",append=TRUE)
      if((((2*m1+1)/2) %in% chr.edges)){
        chr <- chr + 1
      }
    }
  }
  #Annotate the last marker on the last chromosome
  if(!missing(chr.edges)) cat(paste(namez[length(namez)],"\tCHR",chr,"\n",sep=""),file="node_attributes.txt",append=TRUE)
  if(!missing(chr.edges)) cat("Wrote attributes for",length(namez),"markers on",chr,"chromosomes\n")
}

QTLnetwork <- function(CTLobject, lod.threshold = 5, chr.edges, filename){
  idx   <- 1
  edges <- 0
  chr   <- 1
  write.marker.attributes(attr(CTLobject,"qtl"), chr.edges)
  if(missing(filename)){
    filename <- "network_qtl.sif"
    cat("",file=filename)
  }
  namez <- colnames(attr(CTLobject,"qtl"))
  #Connections between genetic markers
  for(m1 in 1:(length(namez)-1)){
    if(missing(chr.edges)){
      cat(paste(namez[m1],"\tCHR\t",namez[m1+1],"\t",1,"\n",sep=""),file=filename,append=TRUE)
    }else{
      if(!(((2*m1+1)/2) %in% chr.edges)){
        cat(paste(namez[m1],"\tCHR",chr,"\t",namez[m1+1],"\t",1,"\n",sep=""),file=filename,append=TRUE)
      }else{
        chr <- chr + 1
      }
    }
  }
  for(x in apply(attr(CTLobject,"qtl"),1,function(x){which(x > lod.threshold)})){
    for(marker in names(x)){
      cat(paste(ctl.name(CTLobject[[idx]]),"\tQTL\t",marker,"\t",round(attr(CTLobject,"qtl")[idx,marker],digits=1),"\n",sep=""),file=filename,append=TRUE)
      edges <- edges+1
    }
    idx <- idx + 1
  }
  cat("Wrote",edges,"edges to",filename,"\n")
}

CTLnetwork <- function(CTLobject, lod.threshold = 5, chr.edges, add.qtl = TRUE, verbose = FALSE){
  cat("",file="network_full.sif")
  cat("",file="network_summary.sif")
  cnt <- 0
  edge_p_count <- NULL
  edge_m_count <- NULL
  for(CTL in CTLobject){
    for(x in 1:ncol(CTL$ctl)){
      if(is.null(CTL$p)) stop("No permutations found, need permutations to determine likelihood\n")
      for(id in which(abs(CTL$l[,x]) > lod.threshold)){
        edge_name <- paste(ctl.name(CTL),"CTL",ctl.names(CTLobject)[x],sep="\t")
        if(edge_name %in% edge_p_count[,1]){
          edgeid <- which(edge_p_count[,1] %in% edge_name)
          if(verbose) cat("Duplicate edge",edgeid,":", edge_name,"\n")
          edge_p_count[edgeid,2] <- as.numeric(edge_p_count[edgeid,2])+abs(CTL$l[id,x])
          edge_p_count[edgeid,3] <- as.numeric(edge_p_count[edgeid,3])+1
        }else{
          edge_p_count <- rbind(edge_p_count,c(edge_name, CTL$l[id,x], 1))
        }
        cat(ctl.name(CTL),"CTL",rownames(CTL$ctl)[id],CTL$l[id,x],colnames(CTL$ctl)[x],"\n",file="network_full.sif",append=TRUE)
        cnt <- cnt+1
      }
    }
    if(verbose)cat("Analysed phenotype:",attr(CTL,"name"),"\n")
  }
  cat("Wrote",cnt,"CTL edges to network_full.sif\n")
  if(!is.null(edge_p_count)){
    
    for(x in 1:nrow(edge_p_count)){
      cat(edge_p_count[x,1],"\t",round(as.numeric(edge_p_count[x,2]),digits=1),"\t",edge_p_count[x,3],"\n",file="network_summary.sif",append=TRUE,sep="")
    }
    cat("Wrote",nrow(edge_p_count),"CTL edges to network_summary.sif\n")
  }
  if(add.qtl) QTLnetwork(CTLobject, lod.threshold, chr.edges, "network_full.sif")
  if(add.qtl) QTLnetwork(CTLobject, lod.threshold, chr.edges, "network_summary.sif")
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
