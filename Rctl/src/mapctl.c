/******************************************************************//**
 * \file Rctl/src/mapctl.c
 * \brief Implementation of functions related to CTLmapping
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "mapctl.h"

#ifdef _OPENMP
  #include <omp.h>
  #define CSTACK_DEFNS 7                                          /* http://stats.blogoverflow.com/2011/08/using-openmp-ized-c-code-with-r/ */

  // extern uintptr_t R_CStackLimit; /* C stack limit */
  // extern uintptr_t R_CStackStart; /* Initial stack address */

  void R_openmp(int* nthr, int* ni, int* ny, double* x, double* ym, double* res) {
    //R_CStackLimit=(uintptr_t) - 1;                              /* http://stats.blogoverflow.com/2011/08/using-openmp-ized-c-code-with-r/ */

    int rthreads = (int)(*nthr);                                  /* Requested number of threads */
    int nitems = (int)(*ni);                                      /* Number of items todo */
    int ylength = (int)(*ny);                                     /* First dimension of y */
    int npthread = ceil((float)nitems / (float)rthreads);         /* Number we do in every thread */
    int tid;                                                      /* private(tid) */
    double** y = asdmatrix(nitems, ylength, ym);                  /* y matrix */

    #pragma omp parallel private(tid) shared(x, y, res) num_threads(rthreads)
    {
      tid = omp_get_thread_num();                                                 /* Obtain thread number */
      int start = npthread * tid;
      int stop = (int)fmin((float)(start + (int)npthread), (float)nitems);
      if (start < stop) {
        //info("I am thread = %d/%d -> [%d,%d]\n", tid, rthreads, start, stop);             /* Echo thread information */
        for(int j = start; j < stop; j++) {                                                 /* Do work we are assigned */
          res[j] = correlation(x, y[j], ylength, false);
        }
      }
    }  /* All threads join master thread and disband */
    free(y);
  }
#else
  void R_openmp(int* nthr, int* ni, int* ny, double* x, double* ym, double* res) {
    info("%s\n", "Unfortunately, openMP is not supported on your platform");
    info("%s\n", "Using a basic for loop on 1 thread");
    int nitems = (int)(*ni);                                      /* Number of items todo */
    int ylength = (int)(*ny);                                     /* First dimension of y */
    double** y = asdmatrix(nitems, ylength, ym);                  /* y matrix */
    for(int j = 0; j < nitems; j++) {                                                 /* Do work we are assigned */
      res[j] = correlation(x, y[j], ylength, false);
    }
    free(y);
  }
#endif

double** mapctl(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                bool doperms, int nperms, int nthreads, bool adjust, bool verbose){

  info("Phenotype %lu: Mapping", (unsigned long)(phenotype+1));
  clvector* genoenc = getGenotypes(genotypes, false);
  size_t i;  
  double** ctls;
  double*  perms;
  double** scores = ctleffects(phenotypes, genotypes, phenotype, genoenc, nthreads, verbose);
  if(!doperms){
    info(", toLOD%s\n", "");  // Exact calculation can be used
    ctls = toLODexact(scores, genoenc, genotypes.nmarkers, phenotypes.nphenotypes, adjust);
  }else{
    info(", Permutation%s", "");
    perms = permute(phenotypes, genotypes, phenotype, genoenc, nperms, nthreads, false);
    info(", toLOD%s\n", "");
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
    err("Cannot scan phenotype %lu out of %lu phenotypes provided", (unsigned long)(phenotype+1), (unsigned long)phenotypes.nphenotypes);
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
          if(verbose) info("Marker %lu, genotype %lu has less then three elements (%lu)\n", (unsigned long)m+1, (unsigned long)g, (unsigned long)idx.nelements);
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
     // warning("Marker %d only has a single genotype\n", m+1);
      dcors[m] = newdvector(phenotypes.nphenotypes);  /*!< Empty Chi^2 values */
    }
  }
  return dcors;
}
