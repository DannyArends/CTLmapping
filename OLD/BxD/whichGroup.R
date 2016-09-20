

whichGroup <- function(x, list){
  return(as.character(list[which(list[,1] == x),2]))
}

#setwd("D:/Github/CTLmapping/examples/BXD/data")         # Where is the data ?
#mlist <- read.table("genes.txt",sep="\t")
