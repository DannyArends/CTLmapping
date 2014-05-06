/******************************************************************//**
 * \file Rctl/src/rmapctl.h
 * \brief Definition of the interfaces to R, and the updateR helper function
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __RMAPCTL_H__
    #define __RMAPCTL_H__

    #include "ctl.h"
    #include "correlation.h"
    #include "permutation.h"    
    #include "sort.h"

    /** Function to 'update' R, checks user input and can flushes console. */
    void     updateR(bool flush);
    /** R interface to perform a CTL scan and permutations on phenotype 'phenotype' */
    void     R_mapctl(int* nind, int* nmar, int* nphe, int* geno, double* pheno,
                      int* p, int *nperms, int* permt, double* dcor, 
                      double* perms, double* res, int* verb);

  #endif //__RMAPCTL_H__
#ifdef __cplusplus
  }
#endif

