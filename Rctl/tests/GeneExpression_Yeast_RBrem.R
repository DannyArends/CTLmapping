#Analysis script for the Yeast GeneExpression
# Paper senior author: RBrem 20??, Published in: Plos
#
# (C) 2012 Danny Arends - GBIC - University of Groningen
#   - Files are distributed with the Rctl package in Rctl/tests
# Start by: 
# setwd("C:/github/CTLmapping/Rctl/tests")
library(ctl)

if(has_rqtl()){
  require(qtl)
  probeannot <- read.csv("GeneExpression_Yeast_Annotation.txt",sep="\t",header=TRUE,row.names=1)
  cross <- read.cross("csvr",file="GeneExpression_Yeast_RBrem.csvr",geno=c("AA","AB"))
  cross <- convert2riself(cross)

  CTL <- CTLscan.cross(cross, 1:5, verbose=T)
  image(CTL)
  plot(CTL)

  nodes <- CTLnetwork(CTL, significance = 0.1)
}else{
  cat("Test skipped because of missing R/qtl package")
}
