#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
      
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
  size_t bufSize = 20 * sizeof(char);
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

  // Read the file content as character
  char* inputAchar = getFilecontent(inputAfn);
  char* inputBchar = getFilecontent(inputBfn);
  char* outputCchar = getFilecontent(outputCfm);

  size_t aol, alc; // a on line, a line count
  size_t bol, blc; // b on line, b line count
  size_t col, clc; // c on line, c line count

  // Convert the file content to double vectors
  double* A = contentToDoubles(inputAchar, &aol, &alc);
  printf("A on line: %zu, A lines: %zu\n", aol, alc);
  double* B = contentToDoubles(inputBchar, &bol, &blc);
  printf("B on line: %zu, B lines: %zu\n", bol, blc);
  double* Co = contentToDoubles(outputCchar, &col, &clc);
  printf("Co on line: %zu, Co lines: %zu\n", col, clc);

  size_t nt = 10;                       // number of times to execute for timing
  size_t nd = 100000000;                // rounding constant
  double accuracy = 0.00001;            // rounding accuracy
  double found, expected;               // rounded down values for comparison
  double sum_time = 0, var_time = 0;    // used in timing the code

  size_t i,j,k, m, mm, n, p, t;
  double sab,sa,sb,saa,sbb;
  m = aol;  // m is the shared dimension between A and B (individuals)
  n = alc;  // n is the dimension unique to A
  p = blc;  // p is the dimension unique to B

  // Allocate the matrix holding the results
  double *C = malloc((n*p) * sizeof(double));
  double *time = malloc((nt) * sizeof(double));

  for (t = 0; t < nt; t++) {  // Run nt times to get a time estimate
    clock_t begin = clock();
    for (i = 0; i < n; i++) {
      for (j = 0; j < p; j++) {
        sa = 0.0;
        sb = 0.0;
        saa = 0.0;
        sbb = 0.0;
        sab = 0.0;

        mm=m;

        for (k = 0; k < m; k++) {
          if ((A[i*m+k] > 2.0) || (B[j*m+k] > 2.0)) {   // Missing value
            mm--;                                       // Decrease the effective number of individuals (mm)
          } else {
            sab += A[i*m+k]*B[j*m+k];
            sa  += A[i*m+k];
            sb  += B[j*m+k];
            saa += A[i*m+k]*A[i*m+k];
            sbb += B[j*m+k]*B[j*m+k];
          }
        }
        C[i*p+j] = (sab-sa*sb/  \
                   (double)mm)/  \
                   (sqrt(saa-sa*sa/  \
                   (double)mm)*sqrt(sbb-sb*sb/  \
                   (double)mm));
      }
    }
    clock_t end = clock();

    // Compare results to the output calculated by R (our gold standard)
    for (i = 0; i < n; i++) {
      for (j = 0; j < p; j++) {
        found = roundf(C[i*p+j] * nd) / nd;
        expected = roundf(Co[i*p+j] * nd) / nd;
        if(fabs(found) < fabs(expected) - accuracy || fabs(found) > fabs(expected) + accuracy) {
          printf("Warning at %zu,%zu, results not equal, found %f, expected %f\n", i, j, found, expected);
        }
      }
    }

    // Update the time information
    time[t] = (double)(end - begin) / CLOCKS_PER_SEC;
    sum_time += time[t];
  }

  // Compute the time standard deviation
  for (t = 0; t < nt; t++) {
    var_time += pow(time[t] - (sum_time / nt), 2);
  }
  printf("Time: %f +/- %f\n", (sum_time / nt), sqrt(var_time / nt));
}

