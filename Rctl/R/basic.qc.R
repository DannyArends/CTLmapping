#
# ctl.print.R
#
# copyright (c) 2010-2011 Danny Arends and Ritsert C. Jansen
# last modified May, 2012
# first written May, 2012
# 
# Basic QC routines used in the examples of CTL analysis
#

### Section: Input data quality control ###

basic.qc <- function(genotypes, phenotypes){
  png("QC_genetic_markers.png",width=2000,height=1000)
  image(1:ncol(genotypes),1:ncol(genotypes),cor(genotypes=="A",use="pair"),
        main="Correlation genetic markers",xlab="Markers",ylab="Markers",col=gray.colors(75))
  abline(h=(get.chredges+.5),col="white")
  abline(v=(get.chredges+.5),col="white"); box()
  dev.off()

  #Create the Two correlation matrices
  cor_metabolites <- cor(metabolites,use='pair')
  rownames(cor_metabolites) <- colnames(metabolites)
  colnames(cor_metabolites) <- colnames(metabolites)
  cor_individuals <- cor(t(metabolites),use='pair')
  rownames(cor_individuals) <- rownames(metabolites)
  colnames(cor_individuals) <- rownames(metabolites)

  png("QC_Metabolites.png",width=2000,height=1000)
  image(x=1:ncol(metabolites),y=1:ncol(metabolites),cor_metabolites ,breaks=seq(-1,1.05,0.05),col=redblue)
  abline(h=(1:ncol(metabolites)+.5),col="white",lwd=0.1)
  abline(v=(1:ncol(metabolites)+.5),col="white",lwd=0.1); box()  
  dev.off()

  png("QC_Individuals.png",width=1000,height=600)
  image(x=1:nrow(metabolites),y=1:nrow(metabolites),cor_individuals ,breaks=seq(-1,1.05,0.05),col=redblue)
  dev.off()

  png("QC_Histogram.png",width=1000,height=600)
  hist(cor_metabolites, xlab="Correlation",main="Correlation: Metabolites (Red), Individuals (Blue)", breaks=seq(-1,1,0.1), col=rgb(1,0,0,0.5),freq=FALSE)
  hist(cor_individuals, breaks=seq(-1,1,0.1), add=TRUE,col=rgb(0,0,1,0.5),freq=FALSE)
  dev.off()

  top.correlated(cor_metabolites)[1:10,]
}

