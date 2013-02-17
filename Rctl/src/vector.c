#include "structs.h"

/* Allocate a new character vector of length dim */
char* newcvector(size_t dim){
  char* v =  calloc(dim, sizeof(char));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

/* Allocate a new double vector of length dim */
double* newdvector(size_t dim){
  double* v = calloc(dim, sizeof(double));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

/* Allocate a new custom length integer vector of length dim */
clvector newclvector(size_t dim){
  clvector clv;
  clv.nelements = 0;
  clv.data = newivector(dim);
  return clv;
}

/* Allocate a new integer vector of length dim */
int* newivector(size_t dim){
  int* v = calloc(dim, sizeof(int));
  if(v==NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  return v;
}

/* Adds a new integer element n to vector v */
int* addtoivector(int* v, size_t dim, int n){
  int* v1 = realloc((void*)v, (dim+1) * sizeof(int));
  if(v1 == NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  v1[dim] = n;
  return v1;
}

/* Adds a new double element n to vector v */
double* addtodvector(double* v, size_t dim, double n){
  double* v1 = realloc((void*)v, (dim+1) * sizeof(double));
  if(v1 == NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  v1[dim] = n;
  return v1;
}

/* Adds a new character element n to vector v */
char* addtocvector(char* v, size_t dim, char n){
  char* v1 = realloc((void*)v, (dim+1));
  if(v1 == NULL) err("Not enough memory for new vector of dimension %d\n",(dim+1));
  v1[dim] = n;
  return v1;
}

/* Print a custom length vector to the output */
void printclvector(const clvector v){
  size_t r;
  info("[");
  for(r = 0; r < v.nelements; r++){ info("%d",v.data[r]); if(r != (dim-1)) info(", "); }
  info("]\n");
}

void printcvector(const char* v, size_t dim){
  size_t r;
  info("[");
  for(r = 0; r < dim; r++){ info("%c",v[r]); if(r != (dim-1)) info(", "); }
  info("]");
}

void printdvector(const double* v, size_t dim){
  size_t r;
  info("[");
  for(r = 0; r < dim; r++){ info("%.2f",v[r]); if(r != (dim-1)) info(", "); }
  info("]");
}

void printivector(const int* v, size_t dim){
  size_t r;
  info("[");
  for(r = 0; r < dim; r++){ info("%d",v[r]); if(r != (dim-1)) info(", "); }
  info("]");
}

clvector which(const int* v, size_t dim, int type){
  size_t i  = 0;
  clvector clv = newclvector(0);
  for(i = 0; i < dim; i++){ 
    if(v[i] == type){
      clv.data = addtoivector(clv.data, clv.nelements, i);
      clv.nelements++;
    }
  }
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

/* Does the clvector contain val */
bool in(const clvector vector, int val){
  size_t i;
  for(i =0; i< vector.nelements; i++){
    if(val == vector.data[i]) return true;
  }
  return false;
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

/* Possibly slow implementation to get the ranks in a double vector 
   Propagates our missing value (-999)
  e.g.: 5      7    8   3   5
  to:   1.5    2    3   1   1.5
*/
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

/* Fisher-Yates random-range generation */
int* randomizeivector(int* idx, size_t max){
  if(max==2) return idx;
  return randomizeivector(swap(idx, (int)(randnum()*(max-2)), max-1),(max-1));
}

