#
# ctl.statistics.R
#
# copyright (c) 2010 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Jan, 2012
# first written Nov, 2011
# 
# R functions to do transform CTL mapping scores to Pvalues and LOD
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#


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
