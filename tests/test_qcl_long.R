library(qcl)
data(ath.metabolites)
qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes, n.perm=100)

jpeg("pxmmatrix.jpg",w=1024,h=768)
image(qcl_result,grid.col="white")
dev.off()

jpeg("pxpmatrix.jpg",w=1024,h=768)
image(qcl_result,grid.col="white",against="phenotypes")
dev.off()

for(x in 1:24){
  jpeg(paste("phenotype",x,".jpg",sep=""),w=1024,h=768)
  plot(qcl_result[[1]])  
  dev.off()
}
