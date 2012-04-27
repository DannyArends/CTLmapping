#
# ctl.network.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Jan, 2012
# first written Oct, 2011
# 
# Network routines for CTL analysis
#

CTLnetwork <- function(CTLscan, ctl.threshold=0.45, verbose = FALSE){
  cat("",file="network.sif")
  cnt <- 0
  edge_count <- NULL
  node_names <- NULL
  for(CTL in CTLscan){
    for(x in 1:ncol(CTL$ctl)){
      if(is.null(CTL$p)) stop("No permutations found, need permutations to determine likelihood\n")
      for(id in which(abs(CTL$ctl[,x]) > getPermuteThresholds(CTL)[1])){
        edge_name <- paste(attr(CTL$ctl,"name"),"CTL",colnames(CTL$ctl)[x],sep="\t")
        if(edge_name %in% edge_count[,1]){
          edgeid <- which(edge_count[,1] %in% edge_name)
          if(verbose) cat("Duplicate edge",edgeid,":", edge_name,"\n")
          edge_count[edgeid,2] <- as.numeric(edge_count[edgeid,2])+abs(CTL$l[id,x])
          edge_count[edgeid,3] <- as.numeric(edge_count[edgeid,3])+1
        }else{
          edge_count <- rbind(edge_count,c(edge_name, CTL$l[id,x], 1))
        }
        node_names <- c(node_names,attr(CTL$ctl,"name"),rownames(CTL$ctl)[id])
        cat(attr(CTL$ctl,"name"),"CTL",rownames(CTL$ctl)[id],CTL$ctl[id,x],colnames(CTL$ctl)[x],"\n",file="network.sif",append=TRUE)
        cnt <- cnt+1
      }
    }
    if(verbose)cat("Analysed phenotype:",attr(CTL,"name"),"\n")
  }
  cat("Wrote",cnt,"edges to network.sif\n")
  if(!is.null(edge_count)){
    cat("",file="edge_summary.sif")
    for(x in 1:nrow(edge_count)){
      cat(edge_count[x,1],"\t",round(as.numeric(edge_count[x,2]),digits=1),"\t",edge_count[x,3],"\n",file="edge_summary.sif",append=TRUE,sep="")
    }
    cat("Wrote",nrow(edge_count),"unique edges to edge_summary.sif\n")
  }
  invisible(unique(node_names))
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
write.nodeAttributeFile <- function(nodenames, spotAnnotation){
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
