/******************************************************************//**
 * \file Rctl/src/permutation.h
 * \brief Definition of functions related to permutations
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __PERMUTATION_H__
    #define __PERMUTATION_H__

    #include "ctl.h"
    #include "mapctl.h"
    #include "sort.h"

    /** Perform permutations using Breitling et al. permutations strategy. */
    double*   permute(const Phenotypes phe, const Genotypes geno, size_t p, clvector* genoenc, 
                      size_t np, int nthreads, bool verbose);
    /** Perform permutations using a row-wise permutations strategy. */
    double**  permuteRW(const Phenotypes phe, const Genotypes geno, size_t p, clvector* genoenc, 
                        size_t np, int nthreads, bool verbose);

    /** Estimate a p-value based on permutations */
    double estimate(double val, double* permutations, size_t nperms);

    /** Converts CTL scores to LOD using exact calculations and bonferonni correction. */
    double** toLODexact(double** scores, clvector* genoenc, size_t nmar, size_t nphe, bool adjust);
    /** Converts CTL scores to LOD using Breitling et al. permutations. */
    double**  toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms);
    /** Converts CTL scores to LOD using row-wise permutations. */
    double**  toLODRW(double** scores, double** permutations, size_t nmar, size_t nphe, size_t nperms);
  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

