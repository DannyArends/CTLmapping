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

#Use the top 10% of permutation scores and fit a GPD uppertail distribution, then use the 
#GPD to obtain P-values for outliers, If the GPD estimates a 0 p-value we reduce our value 
#to get the minimum P-value
extrapolateBeyondRange <- function(permvalues, value = 0.6, top = 20){
  require(POT)
  gpd.threshold <- permvalues[(.80*length(permvalues))]
  mle <- fitgpd(permvalues, gpd.threshold, "mle")
  shape <- mle$param["shape"]
  scale <- mle$scale
  loc <- mle$threshold[1]
  dens <- function(x) dgpd(x, loc, scale, shape)
  warn <- FALSE
  prev.value <- value
  while(as.numeric(dens(value))==0 && value > 0){
    #cat("[FIX] Out of range and p=0:", value," ",as.numeric(dens(value)),"\n")
    warn <- TRUE
    value <- value - 0.0001
  }
  as.numeric(dens(value))
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
