library(qcl)
data(ath.metabolites)
library(qtl)
data(multitrait)
qcl_result <- QCLscan.cross(multitrait, 1:24)
qcl_perms <- QCLpermute(ath.metab$genotypes, ath.metab$phenotypes,1:24,300)

for(x in 1:24){
  jpeg(paste("phenotype",x,".jpg",sep=""),w=1024,h=768)
  op <- par(mfrow=c(2,1))
  plot(QCLscantoScanone(multitrait, qcl_result, qcl_perms, pheno.col=x),scanone(multitrait,pheno.col=x),main="QCL versus QTL")
  plotAsStackedHist(qcl_result, qcl_perms, pheno.col=x, main="Phenotypes contributing to QCL")
  dev.off()
}
