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

double* permute(const Phenotypes phe, const Genotypes geno, size_t p, int ngenotypes,int* genoenc, int a, int b, size_t np, bool verbose){
  size_t perm, ph;
  double* scores = newdvector(np);
  if(a == 1 && b == 1){ return scores; }
  for(perm = 0; perm < np; perm++){
    Genotypes g = permutegenotypes(geno);
    double** ctls  = ctleffects(phe, g, p, ngenotypes, genoenc, a, b, false);
    scores[perm] = fabs(matrixmax(ctls, geno.nmarkers, phe.nphenotypes));
    freematrix((void**)ctls   , geno.nmarkers);
    freematrix((void**)g.data , geno.nmarkers);
    if(verbose) info("Done with permutation %d\n", perm);
    updateR(0);
  }
  qsort(scores, np, sizeof(double), d_cmp);
  return scores;
}

double** permuteRW(const Phenotypes phe, const Genotypes geno, size_t p, int ngenotypes,int* genoenc, int a, int b, size_t np, bool verbose){
  size_t perm,ph;
  double** scores = newdmatrix(phe.nphenotypes, np);
  if(a == 1 && b == 1){ return scores; }
  for(perm = 0; perm < np; perm++){
    Genotypes g = permutegenotypes(geno);
    double** ctls  = ctleffects(phe, g, p, ngenotypes, genoenc, a, b, false);
    double** tctls = transpose(ctls, geno.nmarkers, phe.nphenotypes);
    for(ph = 0; ph <  phe.nphenotypes; ph++){
      scores[ph][perm] = fabs(vectormax(tctls[ph], geno.nmarkers));
    }
    freematrix((void**)ctls   , geno.nmarkers);
    freematrix((void**)tctls  , phe.nphenotypes);
    freematrix((void**)g.data , geno.nmarkers);
    if(verbose) info("Done with permutation %d\n", perm);
    updateR(0);
  }
  for(ph = 0; ph <  phe.nphenotypes; ph++){
    double* ph_s =  scores[ph];
    qsort(ph_s, np, sizeof(double), d_cmp);
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

double** toLODRW(double** scores, double** permutations, size_t nmar, size_t nphe, size_t nperms){
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

