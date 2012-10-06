#include "structs.h"

char* newcvector(size_t dim){
  char* v = (char*) calloc(dim, sizeof(char));
  if(v==NULL){
    printf("Not enough memory for new vector of dimension %d\n",(dim+1));
    exit(-1);
  }
  return v;
}

double* newdvector(size_t dim){
  double* v = (double*) calloc(dim, sizeof(double));
  if(v==NULL){
    printf("Not enough memory for new vector of dimension %d\n",(dim+1));
    exit(-1);
  }
  return v;
}

int* newivector(size_t dim){
  int* v = (int*) calloc(dim, sizeof(int));
  if(v==NULL){
    printf("Not enough memory for new vector of dimension %d\n",(dim+1));
    exit(-1);
  }
  return v;
}

int* addtoivector(const int* v, size_t dim, int n){
  int* v1 = newivector(dim+1);
  size_t i;
  for(i = 0;i < dim;i++){
    v1[i] = v[i];
  }
  v1[dim] = n;
  return v1;
}

double* addtodvector(const double* v, size_t dim, double n){
  double* v1 = newdvector(dim+1);
  size_t i;
  for(i = 0;i < dim;i++){
    v1[i] = v[i];
  }
  v1[dim] = n;
  return v1;
}

char* addtocvector(const char* v, size_t dim, char n){
  char* v1 = newcvector(dim+1);
  size_t i;
  for(i = 0;i < dim;i++){
    v1[i] = v[i];
  }
  v1[dim] = n;
  return v1;
}

void printclvector(clvector v){
  size_t r;
  for(r = 0; r < v.nelements; r++){
    printf("%d, ",v.data[r]);
  }
  printf("\n");
}

void printcvector(const char* v, size_t dim){
  size_t r;
  for(r = 0; r < dim; r++){
    printf("%c, ",v[r]);
  }
  printf("\n");
}

void printdvector(const double* v, size_t dim){
  size_t r;
  for(r = 0; r < dim; r++){
    printf("%f, ",v[r]);
  }
  printf("\n");
}

void printivector(const int* v, size_t dim){
  size_t r;
  for(r = 0; r < dim; r++){
    printf("%d, ",v[r]);
  }
  printf("\n");
}

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

double* get(const double* v, clvector idxs){
  size_t i;
  double* v1 = newdvector(idxs.nelements);
  for(i = 0; i < idxs.nelements; i++){
    v1[i] = v[idxs.data[i]]; 
  }
  return v1;
}
