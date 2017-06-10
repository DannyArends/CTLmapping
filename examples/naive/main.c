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
  
  char*  content = malloc((fsize + 1) * sizeof(char));
  if(content == NULL) printf("Not enough memory for new vector of dimension %zu\n",(fsize));
  
  size_t cnt = 0;
  while(cnt < fsize) {
    content[cnt] = (char)fgetc(file);
    //if(content[cnt] == '\n') printf("%d newline\n", cnt);
    //if(content[cnt] == -1) printf("%d WARN\n", cnt);
    cnt++;
  }
  content[fsize] = -1;
  fclose(file);
  printf("File '%s' loaded: %zu, %zu bytes\n", name, cnt, fsize);
  return content;
}

double* contentToDoubles(char* content, size_t* online, size_t* lines){
  size_t cnt = 0, nCnt = 0, vCnt = 0, lCnt = 0;
  size_t bufSize = 12 * sizeof(char);
  char* nBuffer = malloc(bufSize);
  double* values = malloc(vCnt * sizeof(double));
  while(content[cnt] != -1) {
    //printf("%c %d %d\n", content[cnt], cnt, nCnt);
    if(content[cnt] == '\t'){
      nCnt = 0;
      values = realloc(values, (vCnt+1) * sizeof(double));
      values[vCnt] = atof(nBuffer);
      vCnt++; // numeric value count
      memset(nBuffer, 0, bufSize);
    }else {
      if(content[cnt] == '\n' || content[cnt] == '\0'){
        nCnt = 0;
        values = realloc(values, (vCnt+1) * sizeof(double));
        values[vCnt] = atof(nBuffer);
        vCnt++; // numeric value count
        lCnt++; //line count
        memset(nBuffer, 0, bufSize);
      }else{
        nBuffer[nCnt] = content[cnt];
        nCnt = nCnt + 1; // N char in number buffer count
      }
    }
    cnt = cnt + 1;
  }

  (*lines) = lCnt;          // Number of lines in the file
  (*online) = vCnt / lCnt;  // Number of values on a line
  printf(" %s parsed %zu values from %zu bytes\n", nBuffer, vCnt, cnt);
  return(values);
}

int main(int argc, char **argv){
  char* inputAfn  = "../data/inputA.txt";
  char* inputBfn  = "../data/inputB.txt";
  char* outputCfm = "../data/outputC.txt";

  char* inputAchar = getFilecontent(inputAfn);
  char* inputBchar = getFilecontent(inputBfn);
  char* outputCchar = getFilecontent(outputCfm);
  size_t aol, alc; // a on line, a line count
  size_t bol, blc; // b on line, b line count
  size_t col, clc; // c on line, c line count

  double* A = contentToDoubles(inputAchar, &aol, &alc);
  printf("A on line: %zu, A lines: %zu\n", aol, alc);
  double* B = contentToDoubles(inputBchar, &bol, &blc);
  printf("B on line: %zu, B lines: %zu\n", bol, blc);
  double* Co = contentToDoubles(outputCchar, &col, &clc);
  printf("Co on line: %zu, Co lines: %zu\n", col, clc);

  int i,j,k, m, mm, n, p;
  double sab,sa,sb,saa,sbb;

  m = aol;      // m is the shared dimension between A and B (individuals)
  n = alc;      // n is the dimension unique to A
  p = blc;      // p is the dimension unique to B

  double *C = malloc((n*p) * sizeof(double));
  
  for (i=0; i<n; i++) {
    for (j=0; j<p; j++) {
      sab=0.0;
      sa=0.0;
      sb=0.0;
      saa=0.0;
      sbb=0.0;

      mm=m;

      for (k=0; k<m; k++) {
        if ((A[i*m+k] > 2.0) || (B[j*m+k] > 2.0)) {

          mm--;

        } else {

          sab += A[i*m+k]*B[j*m+k];
          sa  += A[i*m+k];
          sb  += B[j*m+k];
          saa += A[i*m+k]*A[i*m+k];
          sbb += B[j*m+k]*B[j*m+k];

        }
      }
      //printf("%d %d, sab: %f, sa:%f, sb:%f, saa:%f, sbb:%f\n",i, j, sab, sa, sb,saa,sbb);
      C[i*p+j] = (sab-sa*sb/  \
                 (double)mm)/  \
                 (sqrt(saa-sa*sa/  \
                 (double)mm)*sqrt(sbb-sb*sb/  \
                 (double)mm));

    }
  }

  for (i=0; i<n; i++) {
    for (j=0; j<p; j++) {
        if(j < 5) printf("%f == %f\t", C[i*p+j],  Co[i*p+j]);
    }
    printf("\n");
  }

}
