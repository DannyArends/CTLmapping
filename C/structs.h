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

    double**   newdmatrix(size_t rows, size_t cols);
    int**      newimatrix(size_t rows, size_t cols);
    
    double**   addtodmatrix(double** matrix, size_t size, size_t cols, double* n);
    int**      addtoimatrix(int** matrix, size_t size, size_t cols, int* n);

    Phenotypes todmatrix(const char* content);
    Genotypes  toimatrix(const char* content);

    void printdmatrix(double** m, size_t rows, size_t cols);
    void printimatrix(int** m, size_t rows, size_t cols);
    
    void freematrix(void** m, size_t rows);

  #endif //__STRUCTS_H__
#ifdef __cplusplus
  }
#endif
