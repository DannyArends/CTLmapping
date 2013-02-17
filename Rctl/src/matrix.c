#include "matrix.h"

double** newdmatrix(size_t rows, size_t cols){
  size_t i;
  double** m = (double**) calloc(rows, sizeof(double*));
  if(m==NULL) err("Not enough memory for new double matrix [%ix%i]\n", rows, cols);
  for(i = 0; i<rows; i++){ m[i]= newdvector(cols); }
  return m;
}

double** asdmatrix(int rows, int cols, double* data){
  int i;
  double** m = (double**) calloc(rows, sizeof(double*));
  if(m==NULL) err("Not enough memory for new double matrix [%ix%i]\n", rows, cols);
  m[0] = data;
  for(i=1; i< rows; i++){ m[i] = m[i-1] + cols; }
  return m;
}

int** asimatrix(int rows, int cols, int* data){
  int i;
  int** m = (int**) calloc(rows, sizeof(int*));
  if(m==NULL) err("Not enough memory for new integer matrix [%ix%i]\n", rows, cols);
  m[0] = data;
  for(i=1; i< rows; i++){ m[i] = m[i-1] + cols; }
  return m;
}

int** newimatrix(size_t rows, size_t cols){
  size_t x;
  int** m = (int**) calloc(rows, sizeof(int*));
  if(m==NULL) err("Not enough memory for new integer matrix [%ix%i]\n", rows, cols);
  for(x = 0; x<rows; x++){ m[x]= newivector(cols); }
  return m;
}

double** addtodmatrix(double** matrix, size_t size, size_t cols, double* n){
  double** m = realloc(matrix, (size+1) * cols * sizeof(double));
  if(m==NULL) err("Not enough memory for new double matrix [%ix%i]\n", size+1, cols);
  m[size] = n;
  return m;
}

int** addtoimatrix(int** matrix, size_t size, size_t cols, int* n){
  int** m = realloc(matrix, (size+1) * cols * sizeof(int));
  if(m==NULL) err("Not enough memory for new integer matrix [%ix%i]\n", size+1, cols);
  m[size] = n;
  return m;
}

void printdmatrix(double** m, size_t rows, size_t cols){
  size_t r,c;
  for(r = 0; r < rows; r++){
    for(c = 0; c < cols; c++){
      if(c > 0) info("\t");
      info("%f",m[r][c]);
    }
    info("\n");
  }
}

void printimatrix(int** m, size_t rows, size_t cols){
  size_t r,c;
  for(r = 0; r < rows; r++){
    for(c = 0; c < cols; c++){
      if(c > 0) info("\t");
      info("%d",m[r][c]);
    }
    info("\n");
  }
}

/* Get the maximum value in matrix m */
double matrixmax(double** m, size_t rows, size_t cols){
  size_t r,c;
  double max = -DBL_MAX;
  for(r = 0; r < rows; r++){
    for(c = 0; c < cols; c++){
      if(m[r][c] > max) max = m[r][c];
    }
  }
  return max;
}

/* Transpose matrix M */
double** transpose(double** m, size_t rows, size_t cols){
  double** nm = newdmatrix(cols,rows);
  int r, c;
  for(r = 0; r < rows; r++){
    for(c = 0; c < cols; c++){
      nm[c][r] = m[r][c];
    }
  }
  return nm;
}

/* Get the maximum value in vector v */
double vectormax(double* v, size_t dim){
  size_t i;
  double max = -DBL_MAX;
  for(i = 0; i < dim; i++){
    if(v[i] > max) max = v[i];
  }
  return max;
}

void freematrix(void** m, size_t rows){
  size_t i;
  for(i = 0; i < rows; i++){
    free(m[i]);
  }
  free(m);
}

