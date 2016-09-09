library(ctl)
data(ath.metabolites)                 # Arabidopsis Thaliana data set
ctlscan <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, phenocol=1:4)
crs <- getCorrelations(ath.metab$genotypes, ath.metab$phenotypes)
cvs <- chiSQN(crs$correlations, crs$samplesize)
pvs <- chiSQtoP(cvs, (length(crs$samplesize)-1))
q("no")
