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

cross$pheno <- apply(cross$pheno, 2, as.numeric)

qtls  <- scanone(cross, pheno.col = 1:10)   # Scan for QTLS
map   <- qtls[,1:2]                         # Get the genetic map
qtls  <- qtls[,-(1:2)]                      # Get the LOD scores

ctls  <- CTLscan.cross(cross, pheno.col = 1:10, qtls = qtls, verbose=TRUE)
nodes <- ctl.lineplot(ctls, map, significance = 0.1) # Line plot all the phenotypes
nodes

