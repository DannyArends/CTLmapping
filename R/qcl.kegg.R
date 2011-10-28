#
# qcl.kegg.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Oct, 2011
# first written Oct, 2011
# 
# Download gene annotation from KEGG for QCL analysis
#


KEGG.pathway.parser <- function(HTMLdata, verbose = FALSE){
  possible <- regexpr("<nobr><a href=\"/kegg-bin/show_pathway(.+?)</tr>", HTMLdata,useBytes = TRUE)
  mydata <- NULL
  while(possible > 0){
    ppathway <- substr(HTMLdata, possible[[1]], possible[[1]]+attr(possible,"match.length")-11)
    start <- regexpr("<td align=\"left\">",ppathway)[[1]]+17
    ppathway <- substr(ppathway,start,nchar(ppathway))
    if(verbose) cat(" - Possible pathway: ",ppathway,"\n")
    mydata <- rbind(mydata,ppathway)
    HTMLdata <- substr(HTMLdata,possible[[1]]+attr(possible,"match.length")-1,nchar(HTMLdata))
    possible <- regexpr("<nobr><a href=\"/kegg-bin/show_pathway(.+?)</tr>", HTMLdata,useBytes = TRUE)
  }
  invisible(mydata)
}

download.KEGG.annotation <- function(genenames, species = "sce", directory = "tmp_kegg", force.download=FALSE, verbose = FALSE){
  if(.has_rcurl){
    require("RCurl")
    if(missing(genenames)) stop("No gene names  to get annotation for supplied")
    cnt <- 1
    p_cnt <- 1
    if(!file.exists(directory)) dir.create(directory)
    cat("",file="node.kegg.attributes")
    annotations <- vector("list",length(genenames))
    for(gene_name in genenames){
      file_out_name <- paste(directory,"/kegg_",gene_name,".txt",sep="")
      if(file.exists(file_out_name) && !force.download){
        HTMLdata <- paste(scan(file_out_name,what="character",quiet=TRUE),collapse=" ")
        if(verbose) cat("Loaded ",cnt,":",gene_name," from ",file_out_name,"\n",sep="")
      }else{
        url <- paste("http://www.kegg.jp/entry/",species,"/",gene_name,sep="")
        HTMLdata <- getURL(url)
        cat(HTMLdata,file=paste(directory,"/kegg_",gene_name,".txt",sep=""))
        Sys.sleep(runif(1)*2)
        if(verbose) cat("Downloaded ",cnt,":",gene_name," from ",url,"\n",sep="")
      }
      mypathwaydata <- KEGG.pathway.parser(HTMLdata, verbose=verbose)
      if(!is.null(mypathwaydata)){
        cat(gene_name,file="node.kegg.attributes",append=TRUE,sep="")
        pathways <- NULL
        for(pathway in mypathwaydata){
          cat("\t",pathway,file="node.kegg.attributes",append=TRUE,sep="")
          pathways <- c(pathways,pathway)
        }
        annotations[[cnt]] <- pathways
        cat("\n",file="node.kegg.attributes",append=TRUE,sep="")
        p_cnt <- p_cnt+1
      }
      cnt <- cnt+1
    }
    cat("Pathways found for ",p_cnt,"/",cnt," genes\n",sep="")
    invisible(annotations)
  }else{
    stop(.has_rcurl_warnmsg)
  }
}
