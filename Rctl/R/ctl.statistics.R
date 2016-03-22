#
# ctl.statistics.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Nov, 2011
# 
# R functions to do transform CTL mapping scores to Pvalues and LOD
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#

CTLsignificant <- function(CTLobject, significance = 0.05, what = c("names","ids")){
  if(any(class(CTLobject)=="CTLscan")) CTLobject = list(CTLobject)
  all_sign <- NULL
  if(length(what) > 1) what = what[1]
  for(x in 1:length(CTLobject)){ #Get all significant CTLs
    p_above <- which(apply(CTLobject[[x]]$ctl,2,function(y){
    any(y > -log10(significance))}))

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
    colnames(all_sign) <- c("trait","marker","trait","lod","dcor")
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
