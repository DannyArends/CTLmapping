/******************************************************************//**
 * \file Rctl/src/ctlio.c
 * \brief I/O functions for the standalone C application
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "ctlio.h"

/* Write a double matrix to file */
void writeout(double** ctls, size_t phenotype, size_t nmar, size_t nphe){
  size_t p,m;
  FILE* file;
  char* filename = (char*)calloc(25, sizeof(char));
  char* buf      = (char*)calloc(8, sizeof(char));
  strcpy(filename,"pheno");
  sprintf(buf,"%i",phenotype);
  strcat(filename, buf);
  strcat(filename,".o");
  file = fopen(filename,"w+");
  for(p = 0; p < nphe; p++){
    for(m = 0; m < nmar; m++){
      if(m > 0) fprintf(file, "\t");
      fprintf(file, "%f", ctls[m][p]);
    }
    fprintf(file, "\n");
  }
  free(filename);
  free(buf);
  fclose(file);
}

void writesummary(const Phenotypes phenotypes, const Genotypes genotypes, const char* fn, double** ctls, 
                  size_t phenotype, size_t nmar, size_t nphe, clvector* genoenc, double cutoff){
  size_t p, m, i;
  FILE* file;
  if(phenotype == 0){ 
    file = fopen(fn,"w+");
    fprintf(file, "Trait\tMarker\tTrait\tLOD");
    for(i = 0; i < genoenc[0].nelements; i++){ fprintf(file, "\tCor_%d", i); }
    fprintf(file, "\n");
  }else{ file = fopen("summary.txt","a+"); }
  for(p = 0; p < nphe; p++){
    for(m = 0; m < nmar; m++){
      if(ctls[m][p] >= -log10(cutoff)){
        fprintf(file, "%d\t%d\t%d\t%.2f", phenotype, m, p, ctls[m][p]);
        double* cors = getCorrelations(phenotypes, genotypes, phenotype, genoenc[m], m, p, false);
        for(i = 0; i < genoenc[m].nelements; i++){ fprintf(file,"\t%.3f", cors[i]); }
        fprintf(file, "\n");
        free(cors);
      }
    }
  }
  fclose(file);
}

size_t filesize(char* name){
  char   ch;
  size_t cnt  = 0;
  FILE*  file = fopen(name, "r");
  if(file == NULL) err("Error opening file: %s\n", name);
  do{
    ch = fgetc(file);
    cnt++;
  }while(ch != EOF);
  fclose(file);
  return cnt;
}

char* getFilecontent(char* name){
  size_t fsize   = filesize(name);
  char*  content = newcvector((fsize+1));
  FILE*  file    = fopen(name, "r");
  size_t cnt     = 0;
  char   ch;
  if(file == NULL){ err("Error opening file: %s\n", name); }
  do{
    ch = fgetc(file);
    content[cnt] = ch;
    cnt++;
  }while(ch != EOF);
  fclose(file);
  if(content[cnt] != '\n') content = addtocvector(content,cnt,'\n');
  content = addtocvector(content,cnt,'\0');
  info("File '%s' loaded: %d bytes\n", name, fsize);
  return content;
}

