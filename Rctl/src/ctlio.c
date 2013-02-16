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

/* Writes the significant elements to the summary file */
void writesummary(double** ctls, double** scores, size_t phenotype, size_t nmar, size_t nphe, double cutoff){
  size_t p, m;
  FILE* file;
  if(phenotype == 0){ 
    file = fopen("summary.txt","w+");
    fprintf(file, "Trait\tMarker\tTrait\tScore\tLOD\n");
  }else{ file = fopen("summary.txt","a+"); }
  for(p = 0; p < nphe; p++){
    for(m = 0; m < nmar; m++){
      if(ctls[m][p] >= -log10(cutoff)){
        fprintf(file, "%d\t%d\t%d\t%f\t%f\n", phenotype, m, p, scores[m][p], ctls[m][p]);
      }
    }
  }
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

/* Get the whole content of a file as a char[] (Ensures the file is \n terminated) */
char* getfilecontent(char* name){
  size_t fsize   = filesize(name);
  char*  content = newcvector((fsize+1));
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
  content = addtocvector(content,cnt,'\0');
  info("File '%s' loaded: %d bytes\n", name, fsize);
  return content;
}

