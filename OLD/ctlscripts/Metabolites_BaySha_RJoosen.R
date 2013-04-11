#Analysis script for the Arabidopsis Bay x Sha Metabolites
# Paper senior author: RJoosen 2012, Published in: Plant physiology
#
# (C) 2012 Danny Arends - GBIC - University of Groningen
#   - Files are distributed with the Rctl package in Rctl/tests
# Start by: 
# setwd("E:/github/Rpackages/CTLmapping/Rctl/inst/ctlscripts")

metabolites <- read.table("Metabolites_BaySha_RJoosen.txt",sep="\t",row.names=1,header=TRUE)
genotypes   <- read.table("Genotypes_BaySha.txt",row.names=1,header=TRUE)

#Create the genetic map from the first 3 lines og the genotypes file
map_info <- apply(t(genotypes[1:3,]),2,as.numeric)
rownames(map_info) <- colnames(genotypes)

#Genotype matrix by dropping the first 3 lines with map information
genotypes <- genotypes[-c(1:3),]
rownames(genotypes) <- as.numeric(rownames(genotypes))

#Match the genotype matrix to the phenotypes
ids <- match(rownames(metabolites),rownames(genotypes))
genotypes <- genotypes[ids,]

#Print a piece of the data how does it look?
metabolites[1:5,1:10]
genotypes[1:5,1:10]
map_info[1:10,1:3]

#Load the library and scan the data
library(ctl)
basic.qc(genotypes, metabolites, map_info)

#QTLscan because we need to deal with the 4 environments
qtls <- QTLscan(genotypes, metabolites, metabolites[,"Environment"])
ctls <- CTLscan(genotypes, metabolites, geno.enc=c("A","B"), have.qtl=qtls$qtl, n.perm = 1250)

#Create comparison QTL / CTL heatmaps using QTLimage() and image.CTLscan()
png("Comparison_QTL_CTL.png",width=2000,height=1000)
  op <- par(mfrow=c(1,2)); op <- par(cex=1.8);
  QTLimage(ctls); image(ctls);
dev.off()

dir.create("img")

#Plot the individual CTL plot.CTLobject()
for(ctl in ctls){
  png(paste("img/CTL_Metabolites_RJ_",ctl.name(ctl),".png",sep=""),width=2000,height=1000); 
    op <- par(mfrow=c(1,2))
    op <- par(cex=1.8)
    plot(ctl,cex.legend=0.7)
    cnt<-1;add_col <- apply(matrix(0,ncol(ctl$ctl),nrow(ctl$ctl)),2,function(x){x <- rep(cnt,length(x));cnt <<- cnt+1;return(x)})
    image(1:ncol(ctl$ctl), 1:nrow(ctl$ctl),add_col, col=topo.colors(nrow(ctl$ctl),alpha=0.8),ylab="Metabolites",xlab="Markers",main="Phenotype @ Marker x Phenotype matrix")
    image(1:ncol(ctl$ctl), 1:nrow(ctl$ctl),t(ctl$l),col=c(rgb(up,up,up,up)), add=TRUE)
    abline(v=(get.chredges+.5),col="white");box()
    abline(h=(1:nrow(ctl$ctl)+.5),col="white");box()
  dev.off()
}
