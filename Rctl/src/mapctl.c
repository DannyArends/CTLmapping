/******************************************************************//**
 * \file Rctl/src/mapctl.c
 * \brief Implementation of functions related to CTLmapping
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "mapctl.h"
#include <omp.h>

void R_openmp(int* nthr, int* ni, double* res) {
  int rthreads = (int)(*nthr);
  int nitems = (int)(*ni);
  int npthread = ceil((float)nitems / (float)rthreads);
  int nthreads, tid;

  #pragma omp parallel private(nthreads, tid) shared(nitems, npthread, res) num_threads(rthreads)
  {
    tid = omp_get_thread_num();                                                 /* Obtain thread number */
    int start = npthread * tid;
    int stop = (int)fmin((float)(start + (int)npthread), (float)nitems);
    info("I am thread = %d/%d -> [%d,%d]\n", tid, rthreads, start, stop);       /* Echo threadinformation */
    for(int i = start; i < stop; i++) {                                         /* Do work we are assigned */
      res[i] = (double)i;                                                       /* Cores write back in shared memory */
    }
    if (tid == 0) {                                                             /* Only master thread does this */
      nthreads = omp_get_num_threads();
      info("Number of threads = (%d|%d)\n", rthreads, nthreads);
    }
  }  /* All threads join master thread and disband */
}

double** mapctl(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                bool doperms, int nperms, int nthreads, bool verbose){

  info("Phenotype %d: Mapping", (phenotype+1));
  clvector* genoenc = getGenotypes(genotypes, false);
  size_t i;  
  double** ctls;
  double*  perms;
  double** scores = ctleffects(phenotypes, genotypes, phenotype, genoenc, nthreads, verbose);
  if(!doperms){
    info(", toLOD\n", "");  // Exact calculation can be used
    ctls = toLODexact(scores, genoenc, genotypes.nmarkers, phenotypes.nphenotypes);
  }else{
    info(", Permutation", "");
    perms = permute(phenotypes, genotypes, phenotype, genoenc, nperms, nthreads, false);
    info(", toLOD\n", "");
    ctls = toLOD(scores, perms, genotypes.nmarkers, phenotypes.nphenotypes, nperms);
    free(perms);
  }
  for(i = 0; i < genotypes.nmarkers; i++){ free(genoenc[i].data); }
  free(genoenc);
  freematrix((void**)scores, genotypes.nmarkers);
  return ctls;
}

double** ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                    clvector* genoenc, int nthreads, bool verbose){
  size_t g, m, ngenotypes;
  clvector idx;
  double*  P1;
  int* nsamples;
  double** cors, **P2M;
  double** dcors = (double**) calloc(genotypes.nmarkers, sizeof(double*));

  if(phenotype >= phenotypes.nphenotypes){
    err("Cannot scan phenotype %d out of %d phenotypes provided", (phenotype+1), phenotypes.nphenotypes);
  }

  for(m = 0; m < genotypes.nmarkers; m++){
    ngenotypes = genoenc[m].nelements;
    if(ngenotypes > 1) {
      nsamples   = newivector(ngenotypes);
      cors       = calloc(ngenotypes, sizeof(double*));
      for(g = 0; g < ngenotypes; g++){
        idx = which(genotypes.data[m], phenotypes.nindividuals, genoenc[m].data[g]);

        if(idx.nelements > 3) {
          P1            = get(phenotypes.data[phenotype], idx);
          P2M           = getM(phenotypes.data, idx, phenotypes.nphenotypes);
          cors[g]       = cor1toN(P1, P2M, idx.nelements, phenotypes.nphenotypes, nthreads, verbose);
          nsamples[g]   = idx.nelements;
          free(P1);                                         // Clear the indexes and phenotype1 data
          freematrix((void**)P2M, phenotypes.nphenotypes);  // Clear phenotype2M data
        } else {
          if(verbose) info("Marker %d, genotype %d has less then three elements (%d)\n", m+1, g, idx.nelements);
        }

        free(idx.data);
        #ifdef USING_R
          updateR(0);       // annoying function call to not crash R
        #endif //USING_R
      }
      dcors[m] = chiSQN(ngenotypes, cors, phenotype, nsamples, phenotypes.nphenotypes);
      freematrix((void**)cors, ngenotypes);         // Clear correlation and samples data 
      free(nsamples);
    } else {
      warning("Marker %d only has a single genotype\n", m+1);
      dcors[m] = newdvector(phenotypes.nphenotypes);  /*!< Empty Chi^2 values */
    }
  }

  return dcors;
}

