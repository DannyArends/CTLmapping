#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __STRUCTS_H__
    #define __STRUCTS_H__

    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>

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

    typedef struct{
      size_t* data;
      size_t  nelements;
    }IndexVector;

    double*  newdvector(size_t dim);
    int*     newivector(size_t dim);
    char*    newcvector(size_t dim);
    double*  addtodvector(double* vector, size_t size, double n);
    int*     addtoivector(int* vector, size_t size, int n);
    char*    addtocvector(char* vector, size_t size, char n);

    double** newdmatrix(size_t rows, size_t cols);
    int**    newimatrix(size_t rows, size_t cols);

    Phenotypes todmatrix(char* content);
    Genotypes  toimatrix(char* content);

    void printdmatrix(double** m, size_t rows, size_t cols);
    void printimatrix(int** m, size_t rows, size_t cols);

    void freevector(void* v);
    void freematrix(void** m, size_t rows);

    IndexVector which(int* marker, size_t nindividuals, int type);
    double*     get(double* phenotype, IndexVector idx);
  #endif //__STRUCTS_H__
#ifdef __cplusplus
  }
#endif