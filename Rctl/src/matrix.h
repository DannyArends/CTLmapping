/******************************************************************//**
 * \file Rctl/src/matrix.h
 * \brief Definition of functions related to 3D matrices
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __MATRIX_H__
    #define __MATRIX_H__

    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <float.h>
    #include "ctl.h"
    #include "vector.h"

    /** Allocate a new double matrix. */
    double** newdmatrix(size_t rows, size_t cols);
    /** Allocate a new integer matrix. */
    int**    newimatrix(size_t rows, size_t cols);
    
    /** Add a new row to a double matrix. */
    double** addtodmatrix(double** matrix, size_t size, size_t cols, double* n);
    /** Add a new row to an integer matrix. */
    int**    addtoimatrix(int** matrix, size_t size, size_t cols, int* n);

    /** Convert a double vector to a double matrix. */
    double** asdmatrix(int rows, int cols, double* data);
    /** Convert an integer vector to an integer matrix. */
    int**    asimatrix(int rows, int cols, int* data);

    /** Print a double matrix. */
    void printdmatrix(double** m, size_t rows, size_t cols);
    /** Print an integer matrix. */
    void printimatrix(int** m, size_t rows, size_t cols);

    /** Get by index columns from matrix m. */
    double** getM(double** m, clvector idxs, size_t length);

    /** Get the maximum value in matrix m */
    double matrixmax(double** m, size_t rows, size_t cols);

    /** Get the transposition of matrix m. */
    double** transpose(double** m, size_t rows, size_t cols);

    /** Deallocate matrix m. */
    void freematrix(void** m, size_t rows);

  #endif //__MATRIX_H__
#ifdef __cplusplus
  }
#endif
