/******************************************************************//**
 * \file Rctl/src/permutation.c
 * \brief Implementation of functions related to permutations
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
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

double* permute(const Phenotypes phe, const Genotypes geno, size_t p, clvector* genoenc, 
                size_t np, bool verbose){
  size_t perm;
  double* scores = newdvector(np);

  for(perm = 0; perm < np; perm++){
    Genotypes g = permutegenotypes(geno);
    double** ctls  = ctleffects(phe, g, p, genoenc, false);
    scores[perm] = fabs(matrixmax(ctls, geno.nmarkers, phe.nphenotypes));
    freematrix((void**)ctls   , geno.nmarkers);
    freematrix((void**)g.data , geno.nmarkers);
    if(verbose) info("Done with permutation %d\n", perm);
    #ifdef USING_R
      updateR(0);
    #endif //USING_R
  }
  qsort(scores, np, sizeof(double), d_cmp);
  return scores;
}

double** permuteRW(const Phenotypes phe, const Genotypes geno, size_t p, clvector* genoenc, 
                   size_t np, bool verbose){
  size_t perm,ph;
  double** scores = newdmatrix(phe.nphenotypes, np);

  for(perm = 0; perm < np; perm++){
    Genotypes g = permutegenotypes(geno);
    double** ctls  = ctleffects(phe, g, p, genoenc, false);
    double** tctls = transpose(ctls, geno.nmarkers, phe.nphenotypes);
    for(ph = 0; ph <  phe.nphenotypes; ph++){
      scores[ph][perm] = fabs(vectormax(tctls[ph], geno.nmarkers));
    }
    freematrix((void**)ctls   , geno.nmarkers);
    freematrix((void**)tctls  , phe.nphenotypes);
    freematrix((void**)g.data , geno.nmarkers);
    if(verbose) info("Done with permutation %d\n", perm);
    #ifdef USING_R
      updateR(0);
    #endif //USING_R
  }
  for(ph = 0; ph <  phe.nphenotypes; ph++){
    double* ph_s =  scores[ph];
    qsort(ph_s, np, sizeof(double), d_cmp);
    scores[ph] = ph_s;
  }
  return scores;
}

/* Get the index in a permutation row for which permutation[i] >= val. */
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

double** toLODexact(double** scores, clvector* genoenc, size_t nmar, size_t nphe){
  double** ctls = newdmatrix(nmar, nphe);
  size_t p,m;
  double pval;
  for(m = 0; m < nmar; m++){
    size_t Dof = (genoenc[m].nelements-1);
    for(p = 0; p < nphe; p++){
      pval = chiSQtoP(scores[m][p], Dof);
      if(pval > 1 || pval < 0) err("p-value '%.8f' not in range [0, 1]\n", pval);
      pval *= nmar * nphe;
      if(pval >= 1){
        ctls[m][p] = 0.0;
      }else{ ctls[m][p] = fabs(log10(pval)); }
    }
    #ifdef USING_R
      updateR(0);
    #endif //USING_R
  }
  return ctls;
}

double** toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms){
  double** ctls = newdmatrix(nmar, nphe);
  size_t p,m;
  for(m = 0; m < nmar; m++){
    for(p = 0; p < nphe; p++){
      ctls[m][p] = estimate(scores[m][p], permutations, nperms);
    }
    #ifdef USING_R
      updateR(0);
    #endif //USING_R
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
    #ifdef USING_R
      updateR(0);
    #endif //USING_R
  }
  return ctls;
}

