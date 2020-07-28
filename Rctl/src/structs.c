/******************************************************************//**
 * \file Rctl/src/structs.c
 * \brief Implementation of functions related to the Genotype and Phenotype structures
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "structs.h"

clvector* getGenotypes(const Genotypes geno, bool verbose){
  size_t m, i;
  size_t nmar = geno.nmarkers;
  clvector* genoenc = calloc(nmar, sizeof(clvector));
  for(m = 0; m < nmar; m++){
    genoenc[m] = newclvector(0);
    for(i = 0; i < geno.nindividuals; i++){
      if(!CHECKNA(geno.data[m][i]) && !in(genoenc[m], geno.data[m][i])){
        if(verbose) info("Found genotype: %d at marker %d ind %d\n",geno.data[m][i], m, i);
        genoenc[m].data = addtoivector(genoenc[m].data, genoenc[m].nelements, geno.data[m][i]); 
        genoenc[m].nelements++;
      }
    }
  }
  return genoenc;
}

bool checkRow(int ccol, int nrows, int* colcnt){
  if(!(*colcnt)){
    (*colcnt) = ccol;
    return true;
  }else if(ccol != (*colcnt)){
    info("Wrong number of columns on line %d\n",nrows);
    return false;
  }
  return true;
}

Phenotypes toPhenotypes(char* content){
  char*    num    = newcvector(0);
  int      colcnt = 0;
  int      ccol   = 0;
  int      nrows  = 0;
  double** matrix = newdmatrix(0, 0);
  double*  row    = newdvector(0);
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
      if(checkRow(ccol, nrows, &colcnt)){
        num = addtocvector(num, l, '\0');
        row = addtodvector(row, ccol-3, atof(num));
        matrix = addtodmatrix(matrix, nrows, ccol-3, row);
        free(num);
        num  = newcvector(0);
        row  = newdvector(0);
        l    = 0;
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
  if(checkRow(ccol, nrows, &colcnt)){
    num = addtocvector(num, l, '\0');
    row = addtodvector(row, ccol-3, atof(num));
    matrix = addtodmatrix(matrix, nrows, ccol-3, row);
    free(num);
    nrows++;
  }
  info("Parsed into %dx%d matrix\n",nrows,colcnt-2);
  Phenotypes p;
  p.nindividuals = colcnt-2;
  p.nphenotypes = nrows;
  p.data = matrix;
  return p;
}

Genotypes toGenotypes(char* content){
  char*   num = newcvector(0);
  int     colcnt = 0;
  int     ccol   = 0;
  int     nrows  = 0;
  int**   matrix = newimatrix(0, 0);
  int*    row    = newivector(0);
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
      if(checkRow(ccol, nrows, &colcnt)){
        num = addtocvector(num, l, '\0');
        row = addtoivector(row, ccol-3, atoi(num));
        matrix = addtoimatrix(matrix, nrows, ccol-3, row);
        free(num);
        num  = newcvector(0);
        row  = newivector(0);
        l    = 0;
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
  if(checkRow(ccol, nrows, &colcnt)){
    num = addtocvector(num, l, '\0');
    row = addtoivector(row, ccol-3, atoi(num));
    matrix = addtoimatrix(matrix, nrows, ccol-3, row);
    free(num);
    nrows++;
  }
  info("Parsed into %dx%d matrix\n",nrows,colcnt-2);
  Genotypes g;
  g.nindividuals = colcnt-2;
  g.nmarkers = nrows;
  g.data = matrix;
  return g;
}

