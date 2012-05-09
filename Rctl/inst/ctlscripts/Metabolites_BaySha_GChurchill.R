#Analysis script for the Arabidopsis Bay x Sha Metabolites
# Paper senior author: Gary Churchill 2012, Published in: Plos
#
# (C) 2012 Danny Arends - GBIC - University of Groningen
#   - Files are distributed with the Rctl package in Rctl/tests
# Start by: 
# setwd("C:/github/CTLmapping/Rctl/tests")

metabolites <- read.table("Metabolites_BaySha_GChurchill.txt",sep="\t",row.names=1,header=TRUE)
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

#Take only the mean expressions.
#How can we use the StdDev?
metab_means <- grep("Mean",colnames(metabolites))
metabolites <- metabolites[,metab_means]

#Print a piece of the data how does it look?
metabolites[1:5,1:10]
genotypes[1:5,1:10]
map_info[1:10,1:3]

#Load the library and scan the data
library(ctl)
source("Helper_Functions.r")
source("Basic_QC.R",local=TRUE,echo=TRUE)
ctls <- CTLscan(genotypes, metabolites, geno.enc=c("A","B"), n.perm = 1250)

#Create comparison QTL / CTL heatmaps using QTLimage() and image.CTLscan()
png("Comparison_QTL_CTL.png",width=2000,height=1000)
  op <- par(mfrow=c(1,2)); op <- par(cex=1.8);
  QTLimage(ctls); image(ctls);
dev.off()

postscript("QTL_CTL.eps", width = 14.0, height = 6.0)
  op <- par(mfrow=c(1,2)); op <- par(cex=1);
  QTLimage(ctls); image(ctls);
dev.off()


dir.create("img")

#Plot the individual CTL plot.CTLobject()
for(ctl in ctls){
  #postscript(paste("img/CTL_Metabolites_GC_",name(ctl),".eps",sep=""), width = 12.0, height = 6.0,colormodel="rgb")
  png(paste("img/CTL_Metabolites_GC_",name(ctl),".png",sep=""),width=2000,height=1000); 
    op <- par(mfrow=c(1,2))
    op <- par(cex=1.0)
    plot(ctl,cex.legend=0.7)
    abline(v=(get.chredges+.5),col="black");box()
    cnt<-1;add_col <- apply(matrix(0,ncol(ctl$ctl),nrow(ctl$ctl)),2,function(x){x <- rep(cnt,length(x));cnt <<- cnt+1;return(x)})
    image(1:ncol(ctl$ctl), 1:nrow(ctl$ctl),add_col, col=topo.colors(nrow(ctl$ctl),alpha=0.8),ylab="Metabolites",xlab="Markers",main=paste(name(ctl),"@ Marker x Phenotype matrix"))
    image(1:ncol(ctl$ctl), 1:nrow(ctl$ctl),t(ctl$l),col=c(rgb(up,up,up,up)), add=TRUE)
    abline(v=(get.chredges+.5),col="black");box()
    abline(h=(1:nrow(ctl$ctl)+.5),col="white");box()
  dev.off()
}

pxp <- CTLprofiles(ctls,against="phenotypes")

postscript(paste("PXP.eps",sep=""), width = 16.0, height = 14.0,colormodel="rgb")
 png(paste("img/PXP",name(ctl),".png",sep=""),width=2000,height=1000); 
op <- par(cex=1.4)
internal.image(pxp, mainlabel="CTL - Trait x Traits interaction matrix", do.grid=T, grid.col="black")
box()
dev.off()

distancetree <- hclust((dist(pxp)))
plot(as.phylo(distancetree),type="fan")
plot(as.dendrogram(distancetree))

significant <- 1
outmatrix <- NULL
for(mctl in 1:length(ctls)){
  ctl <- ctls[[mctl]]
  for(mrow in 1:nrow(ctl$ctl)){
    if(any(ctl$l[mrow,] > 5)){
      mcol <- which(ctl$l[mrow,] > 5)
      cat(mctl, name(ctl), mrow, rownames(ctl$ctl)[mrow], mcol,"\n",sep='\t')
      outmatrix <- rbind(outmatrix,c(mctl, name(ctl), mrow, rownames(ctl$ctl)[mrow]))
      significant <- significant + 1
    }
  }
}

cat("Numer of traits with CTL: ",length(unique(outmatrix[,1])),"\n")
cat("Numer of CTL: ",length(outmatrix[,1]),"\n")
ctloutmatrix <- NULL
CTLperTrait <- table(as.numeric(outmatrix[,1]))
for(x in 1:length(CTLperTrait)){
  cat(name(ctls[[as.numeric(names(CTLperTrait)[x])]]),CTLperTrait[x],"\n",sep="\t")
  ctloutmatrix <- rbind(ctloutmatrix,c(name(ctls[[as.numeric(names(CTLperTrait)[x])]]),CTLperTrait[x]))
}


