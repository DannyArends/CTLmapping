#include "correlation.h"
#include "structs.h"
 
matrix mapping(matrix phenotypes, matrix genotypes, size_t nind, size_t nmarkers, size_t nphenotypes, size_t phenotype){
  size_t m,p;
  matrix difcormatrix = newmatrix(nmarkers, nphenotypes);
  for(m = 0; m < nmarkers; m++){
    size_t[] ind_aa  = which(genotypes[m],0);
    size_t[] ind_bb  = which(genotypes[m],1);
    double[] pheno_aa = get(phenotypes[phenotype],ind_aa);
    double[] pheno_bb = get(phenotypes[phenotype],ind_bb);
    for(p = 0; p < nphenotypes; p++){
      double cor_aa = correlation(pheno_aa, get(phenotypes[p],ind_aa), nindividuals);
      double cor_bb = correlation(pheno_bb, get(phenotypes[p],ind_bb), nindividuals);
      difcormatrix[m][p] = pow(cor_aa - cor_bb, 2);
    }
  }
  return difcormatrix;
}

