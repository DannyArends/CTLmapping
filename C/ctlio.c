#include "ctlio.h"

size_t filesize(char* name){
  char ch;
  size_t cnt = 0;
  FILE* file = fopen(name, "r");
  if(file == NULL){
    printf("Error opening file: %s\n", name);
    return cnt;
  }
  do{
    ch = fgetc(file);
    cnt++;
  }while(ch != EOF);
  fclose(file);
  return cnt;
}

char* getfilecontent(char* name){
  size_t fsize   = filesize(name);
  char*  content = newcvector(fsize);
  FILE*  file    = fopen(name, "r");
  size_t cnt     = 0;
  char   ch;
  if(file == NULL){
    printf("Error opening file: %s\n", name);
    exit(-1);
  }
  do{
    ch = fgetc(file);
    content[cnt] = ch;
    cnt++;
  }while(ch != EOF);
  fclose(file);
  if(content[cnt] != '\n') content = addtocvector(content,cnt,'\n');
  printf("File '%s' loaded: %d bytes\n", name, fsize);
  return content;
}
