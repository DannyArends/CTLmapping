/******************************************************************//**
 * \file Rctl/src/vector.c
 * \brief Implementation of functions related to 2D vectors
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "structs.h"

char* newcvector(size_t dim){
  char* v =  calloc(dim, sizeof(char));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

double* newdvector(size_t dim){
  double* v = calloc(dim, sizeof(double));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

int* newivector(size_t dim){
  int* v = calloc(dim, sizeof(int));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

clvector newclvector(size_t dim){
  clvector clv;
  clv.nelements = 0;
  clv.data = newivector(dim);
  return clv;
}

char* addtocvector(char* v, size_t dim, char n){
  char* v1 = realloc((void*)v, (dim+1));
  if(v1 == NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  v1[dim] = n;
  return v1;
}

double* addtodvector(double* v, size_t dim, double n){
  double* v1 = realloc((void*)v, (dim+1) * sizeof(double));
  if(v1 == NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  v1[dim] = n;
  return v1;
}

int* addtoivector(int* v, size_t dim, int n){
  int* v1 = realloc((void*)v, (dim+1) * sizeof(int));
  if(v1 == NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  v1[dim] = n;
  return v1;
}

void printcvector(const char* v, size_t dim){
  size_t r;
  info("[", "");
  for(r = 0; r < dim; r++){ info("%c",v[r]); if(r != (dim-1)) info(", ", ""); }
  info("]", "");
}

void printdvector(const double* v, size_t dim){
  size_t r;
  info("[", "");
  for(r = 0; r < dim; r++){ info("%.2f",v[r]); if(r != (dim-1)) info(", ", ""); }
  info("]", "");
}

void printivector(const int* v, size_t dim){
  size_t r;
  info("[", "");
  for(r = 0; r < dim; r++){ info("%d",v[r]); if(r != (dim-1)) info(", ", ""); }
  info("]", "");
}

void printclvector(const clvector v){
  size_t r;
  size_t dim = v.nelements;
  info("[", "");
  for(r = 0; r < dim; r++){ info("%d",v.data[r]); if(r != (dim-1)) info(", ", ""); }
  info("]\n", "");
}

clvector which(const int* v, size_t dim, int e){
  size_t i  = 0;
  clvector clv = newclvector(0);
  for(i = 0; i < dim; i++){ 
    if(v[i] == e){
      clv.data = addtoivector(clv.data, clv.nelements, i);
      clv.nelements++;
    }
  }
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

bool in(const clvector v, int e){
  size_t i;
  for(i =0; i< v.nelements; i++){
    if(e == v.data[i]) return true;
  }
  return false;
}

double vectormax(double* v, size_t dim){
  size_t i;
  double max = -DBL_MAX;
  for(i = 0; i < dim; i++){
    if(v[i] > max) max = v[i];
  }
  return max;
}

/** A random double between 0 and 1 */
double randnum(){
  #ifdef USING_R
    return unif_rand();
  #else
    return rand() / (RAND_MAX+1.0);
  #endif //USING_R
}

/** Swap 2 elements in an integer vectoir */
int* swap(int* idx, int i1, int i2){
  int t = idx[i2];
  idx[i2] = idx[i1];
  idx[i1] = t;
  return idx;
}

double* torank(double* v, size_t dim){
  size_t  i;
  double* r = newdvector(dim);
  size_t  base = 0;
  double  min = DBL_MAX; 
  size_t  e = dim;
  while(e > 0){                // We need to do dim elements;
    clvector ties = newclvector(0);
    for(i = 0; i < dim; i++){  // Scan the array and get the lowest elements
      if(v[i] < min){          // We have a new lowest
        free(ties.data);
        ties = newclvector(0);
        ties.data = addtoivector(ties.data, ties.nelements, i);
        min = v[i];
      }else if(v[i] == min){   // We have an additional tie
        ties.data = addtoivector(ties.data, ties.nelements, i);
      }
    }
    // Ties should now contains the minimum indexes in v
    for(i = 0; i < ties.nelements; i++){
      if(min == MISSING){                      // Propagate missing values (-999)
        r[ties.data[i]] = MISSING;
      }else{
        r[ties.data[i]] = base + (1.0 / ties.nelements);
      }
      v[ties.data[i]] = DBL_MAX;
    }
    e -= ties.nelements;
    base++;
    free(ties.data);
  }
  return r;
}

int* randomizeivector(int* idx, size_t max){
  if(max == 2) return idx;
  return randomizeivector(swap(idx, (int)(randnum()*(max-2)), max-1),(max-1));
}

