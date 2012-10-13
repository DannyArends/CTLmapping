#include "structs.h"

/* Creates the phenotype struct from a string */
Phenotypes tophenotypes(const char* content){
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
        double* t = addtodvector(row, ccol-3, atof(num));
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
        info("Wrong number of columns on line %d\n",nrows);
        rowok = 0;
      }else{
        rowok = 1;
      }
      if(rowok){
        row = addtodvector(row, ccol-3, atof(num));
        num = newcvector(0);
        l = 0;
        matrix = addtodmatrix(matrix, nrows, ccol-3, row);
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
  info("Parsed into %dx%d matrix\n",nrows,colcnt-2);
  Phenotypes p;
  p.nindividuals = colcnt-2;
  p.nphenotypes = nrows;
  p.data = matrix;
  return p;
}

/* Creates the genotype struct from a string */
Genotypes togenotypes(const char* content){
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
        int* t = addtoivector(row, ccol-3, atoi(num));
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
        info("Wrong number of columns on line %d\n",nrows);
        rowok = 0;
      }else{
        rowok = 1;
      }
      if(rowok){
        int* t = addtoivector(row, ccol-3, atoi(num));
        free(row);
        row = t;
        num = newcvector(0);
        l = 0;
        matrix = addtoimatrix(matrix, nrows, ccol-3, row);
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
  info("Parsed into %dx%d matrix\n",nrows,colcnt-2);
  Genotypes g;
  g.nindividuals = colcnt-2;
  g.nmarkers = nrows;
  g.data = matrix;
  return g;
}
