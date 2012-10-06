#include "mapctl.h"

Genotypes permutegenotypes(Genotypes genotypes){
  size_t m,i;
  int* idx = newivector(genotypes.nindividuals);
  for(i = 0; i < genotypes.nindividuals; i++){ idx[i] = i; }

  idx = randomizeivector(idx, genotypes.nindividuals);
  int** newgenodata = newimatrix(genotypes.nmarkers, genotypes.nindividuals);

  for(i = 0; i < genotypes.nindividuals; i++){
    for(m = 0; m < genotypes.nmarkers; m++){
      newgenodata[m][i] = genotypes.data[m][idx[i]];
    }
  }
  free(idx);
  Genotypes g = genotypes;
  g.data = newgenodata;
  return g;
}

double** mapctl(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, int nperms){
  printf("Phenotype %d: Mapping",phenotype);
  double** dcorscores = diffcor(phenotypes, genotypes, phenotype);
  printf(", Permutation");
  double* permutations = permutation(phenotypes, genotypes, phenotype, nperms, 0);
  printf(", toLOD\n");
  double** ctls = toLOD(dcorscores, permutations, genotypes.nmarkers, phenotypes.nphenotypes, nperms);
  freematrix((void**)dcorscores, genotypes.nmarkers);
  free(permutations);
  return ctls;
}

double* permutation(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, size_t nperms, int verbose){
  size_t p;
  double* scores = newdvector(nperms);
  for(p = 0; p < nperms; p++){
    Genotypes g = permutegenotypes(genotypes);
    double** ctls = diffcor(phenotypes, g, phenotype);
    scores[p] = matrixmax(ctls, genotypes.nmarkers, phenotypes.nphenotypes);
    freematrix((void**)ctls   , genotypes.nmarkers);
    freematrix((void**)g.data , genotypes.nmarkers);
    if(verbose) printf("Done with permutation %d\n", p);
  }
  return scores;
}

double getidx(double val, double* permutations, size_t nperms){
  size_t i;
  for(i=0; i < nperms; i++){
    if(permutations[i] > val) return i;
  }
  return (nperms-1);
}

double estimate(double val, double* permutations, size_t nperms){
  return -log10(1.0 - (getidx(val, permutations, nperms)/(double)nperms));
}

double** toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms){
  qsort(permutations, nperms, sizeof(double),d_cmp);
  double** ctls = newdmatrix(nmar, nphe);
  size_t p,m;
  for(m = 0; m < nmar; m++){
    for(p = 0; p < nphe; p++){
      ctls[m][p] = estimate(scores[m][p], permutations, nperms);
    }
  }
  return ctls;
}

double** diffcor(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype){
  size_t m,p,debug = 0;
  double** difcormatrix = newdmatrix(genotypes.nmarkers, phenotypes.nphenotypes);
  for(m = 0; m < genotypes.nmarkers; m++){
    if(debug) printf("Search and allocation marker %d\n", p);
    clvector ind_aa  = which(genotypes.data[m], phenotypes.nindividuals, 0);
    clvector ind_bb  = which(genotypes.data[m], phenotypes.nindividuals, 1);
    double* pheno_aa1 = get(phenotypes.data[phenotype],ind_aa);
    double* pheno_bb1 = get(phenotypes.data[phenotype],ind_bb);
    for(p = 0; p < phenotypes.nphenotypes; p++){
      if(debug) printf("Search and allocation phenotype %d", p);
      double* pheno_aa2 = get(phenotypes.data[p],ind_aa);
      double* pheno_bb2 = get(phenotypes.data[p],ind_bb);
      if(debug) printf(", done");
      double cor_aa = correlation(pheno_aa1, pheno_aa2, ind_aa.nelements);
      double cor_bb = correlation(pheno_bb1, pheno_bb2, ind_bb.nelements);
      if(debug) printf(", correlation done");
      difcormatrix[m][p] = pow(cor_aa - cor_bb, 2);
      if(debug) printf(", cleanup\n");
      free(pheno_aa2);
      free(pheno_bb2);
    }
    if(debug) printf("Marker done\n");
    free(pheno_aa1);
    free(pheno_bb1);
    free(ind_aa.data);
    free(ind_bb.data);
  }
  return difcormatrix;
}
