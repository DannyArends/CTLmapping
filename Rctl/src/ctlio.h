/******************************************************************//**
 * \file Rctl/src/ctlio.h
 * \brief Definition of I/O functions for the standalone C application
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CTLIO_H__
    #define __CTLIO_H__

    #include "ctl.h"
    #include "rmapctl.h"
    #include "correlation.h"
    #include "structs.h"

    void   writeout(double** ctls, size_t phenotype, size_t nmar, size_t nphe);
    /** Writes the significant elements to a summary file */
    void   writesummary(const Phenotypes phenotypes, const Genotypes genotypes, const char* fn, double** ctls, 
                  size_t phenotype, size_t nmar, size_t nphe, clvector* genoenc, double cutoff);

    /** Size in bytes of a file. */
    size_t filesize(char* name);

    /** Get the whole content of a file as a char*.
     *  This function ensures the returned char* is '\n' terminated */
    char*  getFilecontent(char* name);

  #endif //__CTLIO_H__
#ifdef __cplusplus
  }
#endif

