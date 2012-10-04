#ifdef __cplusplus
  extern "C" {
#endif
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
 
  double* newdvector(size_t dim);
  int*    newivector(size_t dim);
  char*   newcvector(size_t dim);
  double* addtodvector(double* vector, size_t size, double n);
  int*    addtoivector(int* vector, size_t size, int n);
  char*   addtocvector(char* vector, size_t size, char n);

  double** newdmatrix(size_t rows, size_t cols);
  int**    newimatrix(size_t rows, size_t cols);
  
  double** todmatrix(char* content);
  int**    toimatrix(char* content);

  void printdmatrix(double** m, size_t rows, size_t cols);
  void printimatrix(int** m, size_t rows, size_t cols);
  
  void freevector(void* v);
  void freematrix(void** m, size_t rows);

  int*    which(int* marker, size_t nindividuals, int type);
  double* get(double* phenotype, size_t nelements, size_t* r);

#ifdef __cplusplus
  }
#endif