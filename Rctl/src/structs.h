#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __STRUCTS_H__
    #define __STRUCTS_H__

    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <float.h>
    #include "vector.h"
    #include "matrix.h"

    typedef struct{
      double** data;
      size_t   nphenotypes;
      size_t   nindividuals;
    }Phenotypes;

    typedef struct{
      int**    data;
      size_t   nmarkers;
      size_t   nindividuals;
    }Genotypes;

    Phenotypes tophenotypes(const char* content);
    Genotypes  togenotypes(const char* content);

  #endif //__STRUCTS_H__
#ifdef __cplusplus
  }
#endif
