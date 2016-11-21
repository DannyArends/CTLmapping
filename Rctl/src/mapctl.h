/******************************************************************//**
 * \file Rctl/src/mapctl.h
 * \brief Definition of functions and includes for CTL mapping
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __MAPCTL_H__
    #define __MAPCTL_H__

    #include "ctl.h"
    #include "rmapctl.h"
    #include "correlation.h"
    #include "permutation.h"
    #include "sort.h"

    /** openMP test function. */
    void R_openmp(int* nthr, int* ni, int* ny, double* x, double* ym, double* res);

    /** Example function that does CTL mapping permutation for a given phenotype. */
    double** mapctl(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                    bool doperms, int nperms, int nthreads, bool verbose);

    /** Get the CTLeffects matrix (Chi square scores) for a given phenotype. */
    double** ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                    clvector* genoenc, int nthreads, bool verbose);

  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

