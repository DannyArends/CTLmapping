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

    double** mapctl(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                    bool doperms, int nperms, bool verbose);

    double** ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                    clvector* genoenc, int alpha, int beta, bool verbose);

  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

