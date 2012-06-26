#CTL with QTL
aa <- genotypes["MSAT5.14"]=="A"
bb <- genotypes["MSAT5.14"]=="B"
cor(metabolites["OHP3.Mean"][aa],metabolites["SOr48.Mean"][aa])
cor(metabolites["OHP3.Mean"][bb],metabolites["SOr48.Mean"][bb])

plot(metabolites["OHP3.Mean"][bb],metabolites["SOr48.Mean"][bb],pch=20,col='red')
points(metabolites["OHP3.Mean"][aa],metabolites["SOr48.Mean"][aa],pch=20,col='blue')

#CTL no QTL
aa <- genotypes["MSAT4.43"]=="A"
bb <- genotypes["MSAT4.43"]=="B"

cor(metabolites["LCRTot.Mean"][aa],metabolites["GSOXLC.Mean"][aa],use="pair")
cor(metabolites["LCRTot.Mean"][bb],metabolites["GSOXLC.Mean"][bb],use="pair")
 
plot(metabolites["LCRTot.Mean"][bb],metabolites["GSOXLC.Mean"][bb],pch=20,col='red')
points(metabolites["LCRTot.Mean"][aa],metabolites["GSOXLC.Mean"][aa],pch=20,col='blue')