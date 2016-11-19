# Analysis script for the Yeast GeneExpression
# Paper senior author: RBrem 20??, Published in: Plos
#
# (C) 2012 Danny Arends - GBIC - University of Groningen
#   - Files are distributed with the Rctl package in Rctl/tests

# Start by: 
# setwd("~/Github/Rpackages/CTLmapping/Rctl/tests")

library(ctl)

probeannot <- read.csv("GeneExpression_Yeast_Annotation.txt", sep="\t", header=TRUE, row.names=1)
cross <- read.cross("csvr", file="GeneExpression_Yeast_RBrem.csvr", geno=c("AA","AB"))
cross <- convert2riself(cross)
cross <- fill.geno(cross)                   # Fill the genotypes
cross <- calc.genoprob(cross)               # Calculate genotype probabilities

ctls  <- CTLscan.cross(cross, phenocol = c(1:5), qtl = FALSE, verbose = TRUE)

for(x in 1:5) {
  qtls <- scanone(cross, pheno.col = x)     # Scan for QTLS
  ctls[[x]]$qtls <- qtls[,3]
}
map   <- qtls[,1:2]                         # Get the genetic map

nodes <- ctl.lineplot(ctls, map, significance = 0.1) # Line plot all the phenotypes
nodes

