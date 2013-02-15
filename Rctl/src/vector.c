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
  int* v1 = newivector(dim+1);
  size_t i;
  for(i = 0;i < dim;i++){
    v1[i] = v[i];
  }
  v1[dim] = n;
  return v1;
}

/* Adds a new double element n to vector v */
double* addtodvector(const double* v, size_t dim, double n){
  double* v1 = newdvector(dim+1);
  size_t i;
  for(i = 0;i < dim;i++){
    v1[i] = v[i];
  }
  v1[dim] = n;
  return v1;
}

/* Adds a new character element n to vector v */
char* addtocvector(const char* v, size_t dim, char n){
  char* v1 = newcvector(dim+1);
  size_t i;
  for(i = 0;i < dim;i++){
    v1[i] = v[i];
  }
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

/* Which integer elements in v are equal to type, returns a clvector of indices */
clvector which(const int* v, size_t dim, int type){
  size_t i;
  size_t length  = 0;
  clvector clv;
  clv.data = newivector(length);
  for(i = 0; i < dim; i++){ 
    if(v[i] == type){
      int* t = addtoivector(clv.data, length, i); 
      free(clv.data);
      clv.data = t;
      length++;
    }
  }  
  clv.nelements = length;
  return clv;
}

/* Get the double elements in v, specified by the indexes in the clvector idxs */
double* get(const double* v, clvector idxs){
  size_t i;
  double* v1 = newdvector(idxs.nelements);
  for(i = 0; i < idxs.nelements; i++){
    v1[i] = v[idxs.data[i]]; 
  }
  return v1;
}

int in(const clvector vector, int val){
  size_t i;
  for(i =0; i< vector.nelements; i++){
    if(val == vector.data[i] && val != -999) return 1;
  }
  return 0;
}

