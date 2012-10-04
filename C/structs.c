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

int* addtoivector(int* vector, size_t size, int n){
  int* v = newivector(size+1);
  size_t i;
  for(i = 0;i < size;i++){
    v[i] = vector[i];
  }
  v[size] = n;
  return v;
}

double* addtodvector(double* vector, size_t size, double n){
  double* v = newdvector(size+1);
  size_t i;
  for(i = 0;i < size;i++){
    v[i] = vector[i];
  }
  v[size] = n;
  return v;
}

char* addtocvector(char* vector, size_t size, char n){
  char* v = newcvector(size+1);
  size_t i;
  for(i = 0;i < size;i++){
    v[i] = vector[i];
  }
  v[size] = n;
  return v;
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
  size_t i;
  int** m = (int**) calloc(rows, sizeof(int*));
  if(m==NULL){
    printf("Not enough memory for new int matrix\n");
    exit(-1);
  }
  for(i = 0; i<rows; i++) {
    m[i]= newivector(cols);
  }
  return m;
}

double** todmatrix(char* content){
  char*   num = newcvector(0);
  int     colcnt = 0;
  int     ccol = 0;
  int     nrows = 0;
  double** matrix = newdmatrix(0, 0);
  double* row = newdvector(0);
  int     rowok = 0;
  int     l = 0;
  do{
    if(content[0] == '\t'){
      row = addtodvector(row, ccol, atof(num));
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
  return matrix;
}

int** toimatrix(char* content){
  char*   num = newcvector(0);
  int     colcnt = 0;
  int     ccol   = 0;
  int     nrows  = 0;
  int** matrix   = newimatrix(0, 0);
  int* row       = newivector(0);
  int     rowok  = 0;
  int     l      = 0;
  do{
    if(content[0] == '\t'){
      row = addtoivector(row, ccol, atoi(num));
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
        row = addtoivector(row, ccol, atoi(num));
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
  return matrix;
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

void freevector(void* v){ free(v); }

void freematrix(void** m, size_t rows){
  size_t i;
  for(i = 0; i<rows; i++) {
    free(m[i]);
  }
  free(m);
}

int* which(int* marker, size_t nindividuals, int type){
  size_t i;
  size_t length  = 0;
  int*   indices = newivector(length);
  for(i = 0; i < nindividuals; i++){ 
    if(marker[i] == type){
      indices = addtoivector(indices, length, i); 
      length++;
    }
  }  
  return indices;
}

double* get(double* phenotype, size_t nelements, size_t* r){
  size_t i;
  double* ph = newdvector(nelements);
  for(i = 0; i < nelements; i++){ ph[i] = phenotype[r[i]]; }
  return ph;
}
