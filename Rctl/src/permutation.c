#include "permutation.h"

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

double* permute(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, int alpha, int gamma, size_t nperms, int verbose){
  size_t p,ph;
  double* scores = newdvector(nperms);
  if(alpha == 1 && gamma == 1){ return scores; }
  for(p = 0; p < nperms; p++){
    Genotypes g = permutegenotypes(genotypes);
    double** ctls  = diffcor(phenotypes, g, phenotype, alpha, gamma);
    scores[p] = fabs(matrixmax(ctls, genotypes.nmarkers, phenotypes.nphenotypes));

    freematrix((void**)ctls   , genotypes.nmarkers);
    freematrix((void**)g.data , genotypes.nmarkers);
    if(verbose) info("Done with permutation %d\n", p);
    updateR(0);
  }
  qsort(scores, nperms, sizeof(double), d_cmp);
  return scores;
}

double** permuteRowWise(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, int alpha, int gamma, size_t nperms, int verbose){
  size_t p,ph;
  double** scores = newdmatrix(phenotypes.nphenotypes, nperms);
  if(alpha == 1 && gamma == 1){ return scores; }
  for(p = 0; p < nperms; p++){
    Genotypes g = permutegenotypes(genotypes);
    double** ctls  = diffcor(phenotypes, g, phenotype, alpha, gamma);
    double** tctls = transpose(ctls, genotypes.nmarkers, phenotypes.nphenotypes);
    for(ph = 0; ph <  phenotypes.nphenotypes; ph++){
      scores[ph][p] = fabs(vectormax(tctls[ph], genotypes.nmarkers));
    }
    freematrix((void**)ctls   , genotypes.nmarkers);
    freematrix((void**)tctls  , phenotypes.nphenotypes);
    freematrix((void**)g.data , genotypes.nmarkers);
    if(verbose) info("Done with permutation %d\n", p);
    updateR(0);
  }
  for(ph = 0; ph <  phenotypes.nphenotypes; ph++){
    double* ph_s =  scores[ph];   
    qsort(ph_s, nperms, sizeof(double), d_cmp);
    scores[ph] = ph_s;
  }
  return scores;
}

double getidx(double val, double* permutations, size_t nperms){
  size_t i;
  for(i=0; i < nperms; i++){
    if(permutations[i] >= fabs(val)-0.00001) return i;
  }
  return (nperms-1);
}

double estimate(double val, double* permutations, size_t nperms){
  return -log10(1.0 - (getidx(val, permutations, nperms)/(double)nperms));
}

double** toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms){
  double** ctls = newdmatrix(nmar, nphe);
  size_t p,m;
  for(m = 0; m < nmar; m++){
    for(p = 0; p < nphe; p++){
      ctls[m][p] = estimate(scores[m][p], permutations, nperms);
    }
    updateR(0);
  }
  return ctls;
}

double** toLODRowWise(double** scores, double** permutations, size_t nmar, size_t nphe, size_t nperms){
  double** ctls = newdmatrix(nmar, nphe);
  size_t p,m;
  for(m = 0; m < nmar; m++){
    for(p = 0; p < nphe; p++){
      ctls[m][p] = estimate(scores[m][p], permutations[p], nperms);
    }
    updateR(0);
  }
  return ctls;
}

