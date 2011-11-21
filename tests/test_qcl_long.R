library(qcl)
data(ath.metabolites)
library(qtl)
data(multitrait)
qcl_result <- QCLscan.cross(multitrait, 1:24)
qcl_perms <- QCLpermute(ath.metab$genotypes, ath.metab$phenotypes,1:24,30)


for(x in 1:24){
  plot(QCLscantoScanone(multitrait, qcl_result, qcl_perms, pheno.col=x),scanone(multitrait,pheno.col=x))
  Sys.sleep(1)
}
