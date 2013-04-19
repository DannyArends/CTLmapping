#
# detect.peaks.R
#
# copyright (c) 2010-2013 - GBIC, Danny Arends and Ritsert C. Jansen
# last modified Apr, 2013
# first written Apr, 2013
# 
# Peak detection algorithm to 'flatten' data above a certain threshold
#

detect.peaks <- function(data, chrEdges = c(1), threshold = 4, verbose = FALSE){
  above   <- as.numeric(data > threshold)
  below   <- -as.numeric((data < (-threshold)))
  pattern <- above + below

  inPeek  <- 0  # Are we in a peak ?
  top     <- 0  # What was our previous TOP score ?
  for(x in 1:length(pattern)){
    if(x %in% chrEdges){  # Crossed the chromosome boundary, reset our 'counters'
      if(inPeek == -1) pattern[topx] <- -2
      if(inPeek == 1) pattern[topx] <- 2
      inPeek <- 0
    }
    if(verbose) cat(x, pattern[x], data[x], inPeek,"\n")
    if(pattern[x] == 1){ # A 1 means above threshold
      if(inPeek == 0){
        inPeek <- 1
        top  <- data[x]
        topx <- x
      }else if(inPeek == 1){ # We're in a peak already check this one is the 'top'
        if(top < data[x]){ top <- data[x] ; topx <- x }
        if(top > data[x]){ pattern[topx] <- 2; top <- data[x]; } # Going down, might be a double peek
      }else{ stop(paste("Positive peak to negative peak at", x)) }
    }else if(pattern[x] == -1){
      if(inPeek == 0){
        inPeek <- (-1)
        top  <- data[x]
        topx <- x
      }else if(inPeek == -1){ # We're in a peak already check this one is the 'top'
        if(top > data[x]){ top <- data[x]; topx <- x } # Its a peak
        if(top < data[x]){ pattern[topx] <- -2; top <- data[x]; } # Going down, might be a double peek
      }else{ stop(paste("Negative peak to positive peak at", x)) }
    }else if(pattern[x] == 0){
      if(inPeek == 1){ inPeek = 0; pattern[topx] <- 2 }
      if(inPeek == -1){ inPeek = 0; pattern[topx] <- -2 }
    }
  }
  pattern
}

# end of detect.peaks.R

