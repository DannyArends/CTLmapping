library(ctl)
data(ath.metabolites)
ctl_result <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, pheno.col=1:3, n.perm=150)

jpeg("pxmmatrix.jpg",w=1024,h=768)
image(ctl_result,grid.col="white")
dev.off()

jpeg("pxpmatrix.jpg",w=1024,h=768)
image(ctl_result,grid.col="white",against="phenotypes")
dev.off()

for(x in 1:3){
  jpeg(paste("phenotype",x,".jpg",sep=""),w=1024,h=768)
  plot(ctl_result[[1]])  
  dev.off()
}
