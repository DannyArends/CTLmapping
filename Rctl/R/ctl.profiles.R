#
# ctl.profiles.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Oct, 2010
# 

#Create the 2 possible CTL matrices: 1) phenotypes versus markers (PxM) and 2) Phenotypes versus phenotypes (PxP)
CTLprofiles <- function(CTLobject, against = c("markers","phenotypes"), significance = 0.05, verbose = FALSE){
  if(missing(CTLobject)) stop("argument 'CTLobject' is missing, with no default")
  mymatrix <- NULL
  mynames <- NULL

  for(p in 1:length(CTLobject)){
    lod <- CTLtoLODvector(CTLobject[[p]], against)
    threshold <- -log10(significance)
    if(max(lod) > threshold){
      mymatrix <- rbind(mymatrix,lod)
      mynames <- c(mynames,ctl.name(CTLobject[[p]]))  
    }
  }
  if(is.null(mymatrix)) { 
    cat("No significant CTLs detected at alpha = ", significance,"\n")
    return(NULL);
  }
  rownames(mymatrix) <- mynames
  if(against[1] == "phenotypes"){
    class(mymatrix) <- c(class(mymatrix),"P2Pmatrix")
  }else{
    class(mymatrix) <- c(class(mymatrix),"P2Mmatrix")
  }
  mymatrix
}

# end of ctl.profiles.R

