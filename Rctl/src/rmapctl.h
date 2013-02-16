#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __RMAPCTL_H__
    #define __RMAPCTL_H__

    #include "ctl.h"
    #include "correlation.h"
    #include "permutation.h"    
    #include "sort.h"
    void     updateR(bool flush);
    void     R_mapctl(int* nind, int* nmar, int* nphe, int* ngeno, int* geno, double* pheno, int* genoenc,
                      int* p, int *nperms, int* a, int* b, int* permt, double* dcor, 
                      double* perms, double* res, int* verb);

  #endif //__RMAPCTL_H__
#ifdef __cplusplus
  }
#endif

