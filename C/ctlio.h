#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CTLIO_H__
    #define __CTLIO_H__
    
    #include "structs.h"
  
    size_t filesize(char* name);
    char*  getfilecontent(char* name);
    
  #endif //__CTLIO_H__
#ifdef __cplusplus
  }
#endif