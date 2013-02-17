#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CTLIO_H__
    #define __CTLIO_H__
    
    #include <string.h>
    #include "ctl.h"
    #include "correlation.h"
    #include "structs.h"
  
    void   writeout(double** ctls, size_t phenotype, size_t nmar, size_t nphe);
/*    void   writesummary(const Phenotypes phenotypes, const Genotypes genotypes, double** ctls, 
                  size_t phenotype, size_t nmar, size_t nphe, size_t ngeno, int* genoenc, double cutoff); */

    size_t filesize(char* name);
    char*  getfilecontent(char* name);
    
  #endif //__CTLIO_H__
#ifdef __cplusplus
  }
#endif

