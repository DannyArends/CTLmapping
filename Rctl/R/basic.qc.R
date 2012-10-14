#
# basic.qc.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written May, 2012
# 
# Basic QC routines used in the examples of CTL analysis
#

### Section: Input data quality control ###

basic.qc <- function(genotypes, phenotypes, map_info){
  png("QC_genetic_markers.png",width=2000,height=1000)
  image(1:ncol(genotypes),1:ncol(genotypes),cor(genotypes=="A",use="pair"),
        main="Correlation genetic markers",xlab="Markers",ylab="Markers",col=gray.colors(75))
  abline(h=(get.chr.edges(map_info)+.5),col="white")
  abline(v=(get.chr.edges(map_info)+.5),col="white"); box()
  dev.off()

  #Create the Two correlation matrices
  cor_phenotypes <- cor(phenotypes,use='pair')
  rownames(cor_phenotypes) <- colnames(phenotypes)
  colnames(cor_phenotypes) <- colnames(phenotypes)
  cor_individuals <- cor(t(phenotypes),use='pair')
  rownames(cor_individuals) <- rownames(phenotypes)
  colnames(cor_individuals) <- rownames(phenotypes)

  png("QC_Phenotypes.png",width=2000,height=1000)
  image(x=1:ncol(phenotypes),y=1:ncol(phenotypes),cor_phenotypes ,breaks=seq(-1,1.05,0.05),col=redblue)
  abline(h=(1:ncol(phenotypes)+.5),col="white",lwd=0.1)
  abline(v=(1:ncol(phenotypes)+.5),col="white",lwd=0.1); box()  
  dev.off()

  png("QC_Individuals.png",width=1000,height=600)
  image(x=1:nrow(phenotypes),y=1:nrow(phenotypes),cor_individuals ,breaks=seq(-1,1.05,0.05),col=redblue)
  dev.off()

  png("QC_Histogram.png",width=1000,height=600)
  hist(cor_phenotypes, xlab="Correlation",main="Correlation: Phenotypes (Red), Individuals (Blue)", 
       breaks=seq(-1,1,0.1), col=rgb(1,0,0,0.5),freq=FALSE)
  hist(cor_individuals, breaks=seq(-1,1,0.1), add=TRUE,col=rgb(0,0,1,0.5),freq=FALSE)
  dev.off()

  invisible(top.correlated(cor_phenotypes))
}

