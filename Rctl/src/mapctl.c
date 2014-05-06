/******************************************************************//**
 * \file Rctl/src/mapctl.c
 * \brief Implementation of functions related to CTLmapping
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "mapctl.h"

double** mapctl(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                bool doperms, int nperms, bool verbose){

  info("Phenotype %d: Mapping", (phenotype+1));
  clvector* genoenc = getGenotypes(genotypes, false);
  size_t i;  
  double** ctls;
  double*  perms;
  double** scores = ctleffects(phenotypes, genotypes, phenotype, genoenc, verbose);
  if(!doperms){
    info(", toLOD\n");  // Exact calculation can be used
    ctls = toLODexact(scores, genoenc, genotypes.nmarkers, phenotypes.nphenotypes);
  }else{
    info(", Permutation");
    perms = permute(phenotypes, genotypes, phenotype, genoenc, nperms, false);
    info(", toLOD\n");
    ctls = toLOD(scores, perms, genotypes.nmarkers, phenotypes.nphenotypes, nperms);
    free(perms);
  }
  for(i = 0; i < genotypes.nmarkers; i++){ free(genoenc[i].data); }
  free(genoenc);
  freematrix((void**)scores, genotypes.nmarkers);
  return ctls;
}

double** ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                    clvector* genoenc, bool verbose){

  size_t g, m, ngenotypes;
  clvector idx;
  double*  P1;
  int* nsamples;
  double** cors, **dcors, **P2M;

  dcors = newdmatrix(genotypes.nmarkers, phenotypes.nphenotypes);
  for(m = 0; m < genotypes.nmarkers; m++){
    ngenotypes = genoenc[m].nelements;
    if(ngenotypes > 1){
      nsamples   = newivector(ngenotypes);
      cors       = calloc(ngenotypes, sizeof(double*));
      for(g = 0; g < ngenotypes; g++){
        idx = which(genotypes.data[m], phenotypes.nindividuals, genoenc[m].data[g]);

        if(idx.nelements > 3){
          P1            = get(phenotypes.data[phenotype], idx);
          P2M           = getM(phenotypes.data, idx, phenotypes.nphenotypes);
          cors[g]       = cor1toN(P1, P2M, idx.nelements, phenotypes.nphenotypes, true);
          nsamples[g]   = idx.nelements;
          free(P1);                                         // Clear the indexes and phenotype1 data
          freematrix((void**)P2M, phenotypes.nphenotypes);  // Clear phenotype2M data
        }

        free(idx.data); 
        updateR(0);
      }
      dcors[m] = chiSQN(ngenotypes, cors, phenotype, nsamples, phenotypes.nphenotypes);
      freematrix((void**)cors, ngenotypes);         // Clear correlation and samples data 
      free(nsamples);
    }
  }
  return dcors;
}

