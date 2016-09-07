library(ctl)
data(ath.metabolites)                 # Arabidopsis Thaliana data set
ctlscan <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, phenocol=1:4)
q("no")
