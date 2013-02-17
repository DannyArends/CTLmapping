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
    
    double*   permute(const Phenotypes phe, const Genotypes geno, size_t p, clvector* genoenc, 
                      int a, int b, size_t np, bool verbose);
    double**  permuteRW(const Phenotypes phe, const Genotypes geno, size_t p, clvector* genoenc, 
                      int a, int b, size_t np, bool verbose);

    double**  toLODexact(double** scores, clvector* genoenc, size_t nmar, size_t nphe);
    double**  toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms);
    double**  toLODRW(double** scores, double** permutations, size_t nmar, size_t nphe, size_t nperms);
  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

