#
# ctl.statistics.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Apr, 2016
# first written Nov, 2011
# 
# R functions to do transform CTL mapping scores to Pvalues and LOD and collect significant results per region
#

CTLregions <- function(CTLobject, mapinfo, phenocol = 1, significance = 0.05, verbose = TRUE) {
  if(any(class(CTLobject)=="CTLscan")) CTLobject = list(CTLobject)
  if(missing(mapinfo)) stop("You need to provide a map object with 'chr' and 'pos' columns")
  if(!all(colnames(mapinfo) %in% c("chr", "pos"))) stop("You need to provide a map object with 'chr' and 'pos' columns")
  p_above <- which(apply(CTLobject[[phenocol]]$ctl, 2, function(y){
    any(y > -log10(significance))
  }))
  if(length(p_above) == 0) stop("No significant CTLs found at threshold = ", significance)
  significant <- ctlscan[[1]]$ctl[, p_above]
  map <- mapinfo[rownames(significant),]
  if(!all(rownames(significant) %in% rownames(map))) stop("CTL markers do not match the provided map")
  regions <- NULL
  for(x in 1:ncol(significant)){
    peeks <- detect.peaks(significant[,x], threshold = -log10(significance))
    for(peek in which(peeks == 2)){
      left <- peek - 1
      right <- peek + 1
      chr <- map[peek, "chr"]
      while(left >= 1 && peeks[left] >= 1 && map[left, "chr"] == chr){
        left <- left - 1
      }
      while(right <= length(peeks) && peeks[right] >= 1 && map[right, "chr"] == chr){
        right <- right + 1
      }
      if(map[left, "chr"] != chr){ pos_s <- 0; }else{ pos_s <- map[left, "pos"]; }
      if(map[right, "chr"] != chr){ pos_e <- map[(right-1), "pos"]; }else{ pos_e <- map[right, "pos"]; }
      if(verbose) {
        cat(ctl.names(CTLobject)[phenocol], "with", colnames(significant)[x], "from marker", left, "to", right, paste0("chr ", chr,":", pos_s,"-", pos_e), "\n")
      }
      regions <- rbind(regions, c(ctl.names(CTLobject)[phenocol], colnames(significant)[x], chr, pos_s, pos_e))
    }
  }
  colnames(regions) <- c("pheno1", "pheno2", "chr", "start", "end")
  return(invisible(data.frame(regions)))
}

CTLsignificant <- function(CTLobject, significance = 0.05, what = c("names", "ids")){
  if(any(class(CTLobject)=="CTLscan")) CTLobject = list(CTLobject)
  all_sign <- NULL
  if(length(what) > 1) what = what[1]
  for(x in 1:length(CTLobject)){ #Get all significant CTLs
    p_above <- which(apply(CTLobject[[x]]$ctl, 2, function(y){
      any(y > -log10(significance))
    }))

    pnames <- colnames(CTLobject[[x]]$ctl)
    mnames <- rownames(CTLobject[[x]]$ctl)

    if(what != "ids"){ p_above <- pnames[p_above] }
    if(length(p_above) > 0){
      for(p in p_above){
        m_above <- which(CTLobject[[x]]$ctl[,p] > -log10(significance))
        if(what != "ids"){ m_above <- mnames[m_above] }
        for(m in m_above){
          if(what[1] == "ids"){
            all_sign <- rbind(all_sign, c(x, m, p, CTLobject[[x]]$ctl[m, p], CTLobject[[x]]$dcor[m, p]) )
          }else{
            all_sign <- rbind(all_sign, c(ctl.name(CTLobject[[x]]), m, p, CTLobject[[x]]$ctl[m, p], CTLobject[[x]]$dcor[m, p]) )
          }
        }
      }
    }
  }
  items <- 0
  if(!is.null(all_sign)){
    all_sign <- as.data.frame(all_sign)
    all_sign[,4] <- round(as.numeric(as.character(all_sign[,4])), digits=2)
    colnames(all_sign) <- c("trait", "marker", "trait", "lod", "dcor")
    items <- nrow(all_sign)
  }
  cat("Found",items,"significant CTLs\n")
  return(all_sign)
}

CTLtoLODvector <- function(CTLscan, against = c("markers","phenotypes")){
  if(missing(CTLscan)) stop("argument 'CTLscan' is missing, with no default")
  if(against[1]=="markers")return(apply(CTLscan$ctl,1,sum))
  if(against[1]=="phenotypes")return(apply(CTLscan$ctl,2,max))
}

CTLscantoScanone <- function(cross, CTLscan){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(missing(CTLscan)) stop("argument 'CTLscan' is missing, with no default")
  
  scores <- CTLtoLODvector(CTLscan)
  scores[which(!is.finite(scores))] <- NA
  lodscorestoscanone(cross, scores)
}

# end of ctl.statistics.R
