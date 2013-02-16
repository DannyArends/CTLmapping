#include "structs.h"

clvector getGenotypes(const Genotypes g){
  size_t m,i;
  clvector clv = newclvector(0);
  for(m = 0; m < g.nmarkers; m++){
    for(i = 0; i < g.nindividuals; i++){
      if(!in(clv, g.data[m][i])){
        // info("Found genotype: %d\n",g.data[m][i]);        
        clv.data = addtoivector(clv.data, clv.nelements, g.data[m][i]); 
        clv.nelements++;
      }
    }
  }
  return clv;
}

/* Creates the phenotype struct from a string */
Phenotypes tophenotypes(char* content){
  char*    num    = newcvector(0);
  int      colcnt = 0;
  int      ccol   = 0;
  int      nrows  = 0;
  double** matrix = newdmatrix(0, 0);
  double*  row    = newdvector(0);
  int      rowok  = 0;
  int      l      = 0;
  do{//    info("%d %d %c\n", ccol, l, content[0]);
    if(content[0] == '\t'){
      if(ccol > 2){
        num = addtocvector(num, l, '\0');
        row = addtodvector(row, ccol-3, atof(num));
      }
      free(num);
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
        num = addtocvector(num, l, '\0');
        row = addtodvector(row, ccol-3, atof(num));
        free(num);
        num = newcvector(0);
        l = 0;
        matrix = addtodmatrix(matrix, nrows, ccol-3, row);
        row = newdvector(0);
        ccol = 0;
        nrows++;
      }
    }
    if(content[0] != ' '){
      num = addtocvector(num, l, content[0]);
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
Genotypes togenotypes(char* content){
  char*   num = newcvector(0);
  int     colcnt = 0;
  int     ccol   = 0;
  int     nrows  = 0;
  int**   matrix = newimatrix(0, 0);
  int*    row    = newivector(0);
  int     rowok  = 0;
  int     l      = 0;
  do{//    info("%d %d %c\n", ccol, l, content[0]);
    if(content[0] == '\t'){
      if(ccol > 2){
        num = addtocvector(num, l, '\0');
        row = addtoivector(row, ccol-3, atoi(num));
      }
      free(num);
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
        num = addtocvector(num, l, '\0');
        row = addtoivector(row, ccol-3, atoi(num));
        free(num);
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

