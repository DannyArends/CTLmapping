#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
      
char* getFilecontent(char* name) {
  FILE*  file = fopen(name, "r");
  if(file == NULL){ printf("Error opening file: %s\n", name); }

  fseek(file, 0L, SEEK_END);
  size_t fsize = ftell(file);
  fseek(file, 0L, SEEK_SET);
  
  char*  content = malloc(fsize * sizeof(char));
  if(content == NULL) printf("Not enough memory for new vector of dimension %d\n",(fsize));
  
  size_t cnt = 0;
  while(cnt < fsize) {
    content[cnt] = (char)fgetc(file);
    if(content[cnt] == '\n') printf("%d newline\n", cnt);
    if(content[cnt] == -1) printf("%d WARN\n", cnt);
    cnt++;
  }
  fclose(file);
  printf("File '%s' loaded: %d, %d bytes\n", name, cnt, fsize);
  return content;
}

double* contentToDoubles(char* content){
  size_t cnt = 0, nCnt = 0, vCnt = 0;
  char* nBuffer = malloc(12 * sizeof(char));
  double* values = malloc(0 * sizeof(double));
  while(content[cnt] != EOF) {
    printf("%c %d %d\n", content[cnt], cnt, nCnt);
    if(content[cnt] == '\t'){
      nCnt = 0;
      vCnt ++;
    }else {
      if(content[cnt] == '\n' || content[cnt] == '\0'){
        nCnt = 0;
        vCnt ++;
      }else{
        nBuffer[nCnt] = content[cnt];
        nCnt = nCnt + 1;
      }
    }
    cnt = cnt + 1;
  }
  printf(" %s parsed %d values from %d bytes\n", nBuffer, vCnt, cnt);
  return(values);
}

int main(int argc, char **argv){
  char* inputAfn  = "../data/inputAb.txt";
  char* inputBfn  = "../data/inputB.txt";
  char* outputCfm = "../data/outputC.txt";

  char* inputAchar = getFilecontent(inputAfn);
  //char* inputBchar = getFilecontent(inputBfn);
  //char* outputCchar = getFilecontent(outputCfm);
  contentToDoubles(inputAchar);
  int i,j,k, m, mm, n, p;
  int sab,sa,sb,saa,sbb;

  double *A, *B, *C;
  
  for (i=0; i<n; i++) {
    for (j=0; j<p; j++) {
      sab=0.0;
      sa=0.0;
      sb=0.0;
      saa=0.0;
      sbb=0.0;

      mm=m;

      for (k=0; k<m; k++) {
        if ((A[i*m+k] > 2.0) || (B[k*p+j] > 2.0)) {

          mm--;

        } else {

          sab+=A[i*m+k]*B[k*p+j];
          sa+=A[i*m+k];
          sb+=B[k*p+j];
          saa+=A[i*m+k]*A[i*m+k];
          sbb+=B[k*p+j]*B[k*p+j];

        }
      }

      C[i*p+j] = (sab-sa*sb/  \
                 (double)mm)/  \
                 (sqrt(saa-sa*sa/  \
                 (double)mm)*sqrt(sbb-sb*sb/  \
                 (double)mm));

    }
  }
}