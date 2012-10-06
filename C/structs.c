#include "structs.h"

double** newdmatrix(size_t rows, size_t cols){
  size_t i;
  double** m = (double**) calloc(rows, sizeof(double*));
  if(m==NULL){
    printf("Not enough memory for new double matrix\n");
    exit(-1);
  }
  for(i = 0; i<rows; i++) {
    m[i]= newdvector(cols);
  }
  return m;
}

int** newimatrix(size_t rows, size_t cols){
  size_t x;
  int** m = (int**) calloc(rows, sizeof(int*));
  if(m==NULL){
    printf("Not enough memory for new int matrix\n");
    exit(-1);
  }
  for(x = 0; x<rows; x++) {
    m[x]= newivector(cols);
  }
  return m;
}

double** addtodmatrix(double** matrix, size_t size, size_t cols, double* n){
  double** m = newdmatrix(size+1, cols);
  size_t i;
  for(i = 0;i < size;i++){
    m[i] = matrix[i];
  }
  m[size] = n;
  return m;
}

int** addtoimatrix(int** matrix, size_t size, size_t cols, int* n){
  int** m = newimatrix(size+1, cols);
  size_t i;
  for(i = 0;i < size;i++){
    m[i] = matrix[i];
  }
  m[size] = n;
  return m;
}

Phenotypes todmatrix(const char* content){
  char*    num    = newcvector(0);
  int      colcnt = 0;
  int      ccol   = 0;
  int      nrows  = 0;
  double** matrix = newdmatrix(0, 0);
  double*  row    = newdvector(0);
  int      rowok  = 0;
  int      l      = 0;
  do{
    if(content[0] == '\t'){
      if(ccol > 2){
        double* t = addtodvector(row, ccol, atof(num));
        free(row);
        row = t;
      }
      num = newcvector(0);
      l = 0;
      ccol++;
    }
    if(content[0] == '\n' || content[0] == '\0'){
      if(!colcnt){
        colcnt = ccol;
        rowok = 1;
      }else if(ccol != colcnt){
        printf("Wrong number of columns on line %d\n",nrows);
        rowok = 0;
      }else{
        rowok = 1;
      }
      if(rowok){
        row = addtodvector(row, ccol, atof(num));
        num = newcvector(0);
        l = 0;
        matrix = addtodmatrix(matrix, nrows, ccol, row);
        row = newdvector(0);
        ccol = 0;
        nrows++;
      }
    }
    if(content[0] != ' '){
      num = addtocvector(num,l,content[0]);
      l++;
    }
    content++;
  }while(content[0] != '\0');
  printf("Parsed into %dx%d matrix\n",nrows,colcnt);
  Phenotypes p;
  p.nindividuals = colcnt-2;
  p.nphenotypes = nrows;
  p.data = matrix;
  return p;
}

Genotypes toimatrix(const char* content){
  char*   num = newcvector(0);
  int     colcnt = 0;
  int     ccol   = 0;
  int     nrows  = 0;
  int**   matrix = newimatrix(0, 0);
  int*    row    = newivector(0);
  int     rowok  = 0;
  int     l      = 0;
  do{
    if(content[0] == '\t'){
      if(ccol > 2){
        int* t = addtoivector(row, ccol, atoi(num));
        free(row);
        row = t;
      }
      num = newcvector(0);
      l = 0;
      ccol++;
    }
    if(content[0] == '\n' || content[0] == '\0'){
      if(!colcnt){
        colcnt = ccol;
        rowok = 1;
      }else if(ccol != colcnt){
        printf("Wrong number of columns on line %d\n",nrows);
        rowok = 0;
      }else{
        rowok = 1;
      }
      if(rowok){
        int* t = addtoivector(row, ccol, atoi(num));
        free(row);
        row = t;
        num = newcvector(0);
        l = 0;
        matrix = addtoimatrix(matrix, nrows, ccol, row);
        row = newivector(0);
        ccol = 0;
        nrows++;
      }
    }
    if(content[0] != ' '){
      num = addtocvector(num,l,content[0]);
      l++;
    }
    content++;
  }while(content[0] != '\0');
  printf("Parsed into %dx%d matrix\n",nrows,colcnt);
  Genotypes g;
  g.nindividuals = colcnt-2;
  g.nmarkers = nrows;
  g.data = matrix;
  return g;
}

void printdmatrix(double** m, size_t rows, size_t cols){
  size_t r,c;
  for(r = 0; r < rows; r++){
    for(c = 0; c < cols; c++){
      if(c > 0) printf("\t");
      printf("%f",m[r][c]);
    }
    printf("\n");
  }
}

void printimatrix(int** m, size_t rows, size_t cols){
  size_t r,c;
  for(r = 0; r < rows; r++){
    for(c = 0; c < cols; c++){
      if(c > 0) printf("\t");
      printf("%d",m[r][c]);
    }
    printf("\n");
  }
}

void freematrix(void** m, size_t rows){
  size_t i;
  for(i = 0; i < rows; i++){
    free(m[i]);
  }
  free(m);
}
