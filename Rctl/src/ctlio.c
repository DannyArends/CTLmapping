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
  file = fopen(filename,"a+");
  for(p = 0; p < nphe; p++){
    for(m = 0; m < nmar; m++){
      fprintf(file, "%f", ctls[m][p]);
      if(m > 0 && m < (nmar-1)) fprintf(file, "\t");
    }
    fprintf(file, "\n");
  }
  free(filename);
  free(buf);
  fclose(file);
}

/* File size of a file */
size_t filesize(char* name){
  char ch;
  size_t cnt = 0;
  FILE* file = fopen(name, "r");
  if(file == NULL){
    info("Error opening file: %s\n", name);
    return cnt;
  }
  do{
    ch = fgetc(file);
    cnt++;
  }while(ch != EOF);
  fclose(file);
  return cnt;
}

/* Get the whole string content of a file */
char* getfilecontent(char* name){
  size_t fsize   = filesize(name);
  char*  content = newcvector(fsize);
  FILE*  file    = fopen(name, "r");
  size_t cnt     = 0;
  char   ch;
  if(file == NULL){
    err("Error opening file: %s\n", name);
  }
  do{
    ch = fgetc(file);
    content[cnt] = ch;
    cnt++;
  }while(ch != EOF);
  fclose(file);
  if(content[cnt] != '\n') content = addtocvector(content,cnt,'\n');
  info("File '%s' loaded: %d bytes\n", name, fsize);
  return content;
}

