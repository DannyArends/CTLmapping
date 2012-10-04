#include "mapctl.h"
 
double** ctlmapping(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype){
  size_t m,p;
  double** difcormatrix = newdmatrix(genotypes.nmarkers, phenotypes.nphenotypes);
  for(m = 0; m < genotypes.nmarkers; m++){
    IndexVector ind_aa  = which(genotypes.data[m], phenotypes.nindividuals,0);
    IndexVector ind_bb  = which(genotypes.data[m], phenotypes.nindividuals,1);
    double* pheno_aa = get(phenotypes.data[phenotype],ind_aa);
    double* pheno_bb = get(phenotypes.data[phenotype],ind_bb);
    for(p = 0; p < phenotypes.nphenotypes; p++){
      double cor_aa = correlation(pheno_aa, get(phenotypes.data[p],ind_aa), phenotypes.nindividuals);
      double cor_bb = correlation(pheno_bb, get(phenotypes.data[p],ind_bb), phenotypes.nindividuals);
      difcormatrix[m][p] = pow(cor_aa - cor_bb, 2);
    }
  }
  return difcormatrix;
}

