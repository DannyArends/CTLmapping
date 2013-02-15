#include "structs.h"

/* Allocate a new character vector of length dim */
char* newcvector(size_t dim){
  char* v = (char*) calloc(dim, sizeof(char));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

/* Allocate a new double vector of length dim */
double* newdvector(size_t dim){
  double* v = (double*) calloc(dim, sizeof(double));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

/* Allocate a new integer vector of length dim */
int* newivector(size_t dim){
  int* v = (int*) calloc(dim, sizeof(int));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

/* Adds a new integer element n to vector v */
int* addtoivector(const int* v, size_t dim, int n){
  int* v1 = (int*)realloc((void*)v, (dim+1) * sizeof(int));
  v1[dim] = n;
  return v1;
}

/* Adds a new double element n to vector v */
double* addtodvector(const double* v, size_t dim, double n){
  double* v1 = (double*)realloc((void*)v, (dim+1) * sizeof(double));
  v1[dim] = n;
  return v1;
}

/* Adds a new character element n to vector v */
char* addtocvector(const char* v, size_t dim, char n){
  char* v1 = (char*)realloc((void*)v, (dim+1) * sizeof(char));
  v1[dim] = n;
  return v1;
}

/* Print a custom length vector to the output */
void printclvector(clvector v){
  size_t r;
  for(r = 0; r < v.nelements; r++){
    info("%d, ",v.data[r]);
  }
  info("\n");
}

void printcvector(const char* v, size_t dim){
  size_t r;
  for(r = 0; r < dim; r++){
    info("%c, ",v[r]);
  }
  info("\n");
}

void printdvector(const double* v, size_t dim){
  size_t r;
  for(r = 0; r < dim; r++){
    info("%f, ",v[r]);
  }
  info("\n");
}

void printivector(const int* v, size_t dim){
  size_t r;
  for(r = 0; r < dim; r++){
    info("%d, ",v[r]);
  }
  info("\n");
}

/* A random double between 0 and 1 */
double randnum(){ return rand() / (RAND_MAX+1.0);}

/* Swap 2 elements in an integer vectoir */
int* swap(int* idx, int i1, int i2){
  int t = idx[i2];
  idx[i2] = idx[i1];
  idx[i1] = t;
  return idx;
}

/* Fisher-Yates random-range generation */
int* randomizeivector(int* idx, size_t max){
  if(max==2) return idx;
  return randomizeivector(swap(idx, (int)(randnum()*(max-2)), max-1),(max-1));
}

