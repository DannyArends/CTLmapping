/*******************************************************************************
* Copyright 2003-2016 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
!  Content:
!    Calculation of covariance and correlation matrices  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM      3       /* Task dimension */
#define N        10000   /* Number of observations */

#define P_THRESHOLD  0.01

/* Exact covariance matrix */
double C[DIM][DIM] = {
    { 16.0,  8.0,  4.0 },
    {  8.0, 13.0, 17.0 },
    {  4.0, 17.0, 62.0 }
};

/* Exact mean vector */
double a[DIM] = { 3.0, 5.0, 2.0 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT x_storage;
    MKL_INT cov_storage;
    MKL_INT cor_storage;
    double x[N][DIM];
    double cov[DIM][DIM], cor[DIM][DIM];
    double mean[DIM];
    double aTmp, rTmp;
    int i, j, errcode;
    int errnums = 0;

    double pval_mean[DIM];
    double pval_cov[DIM][DIM];

    /***** Initializing parameters for Summary Statistics task *****/
    dim         = DIM;
    n           = N;
    x_storage   = VSL_SS_MATRIX_STORAGE_COLS;
    cov_storage = VSL_SS_MATRIX_STORAGE_FULL;
    cor_storage = VSL_SS_MATRIX_STORAGE_FULL;

    for(i = 0; i < dim; i++)
    {
        mean[i] = 0.0;
        for(j = 0; j < dim; j++)
        {
            cov[i][j] = 0;
            cor[i][j] = 0;
        }
    }

    /***** Generate data set using VSL GaussianMV RNG *****/
    errcode = dGenerateGaussianMVData( (double*)x, dim, n, (double*)a, (double*)C );
    CheckVslError(errcode);

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, &n, &x_storage, (double*)x, 0, 0 );
    CheckVslError(errcode);

    /***** Initialization of the task parameters using FULL_STORAGE
        for covariance/correlation matrices *****/
    errcode = vsldSSEditCovCor( task, mean, (double*)cov, &cov_storage,
        (double*)cor, &cor_storage );
    CheckVslError(errcode);

    /***** Compute covariance/correlation matrices using FAST method  *****/
    errcode = vsldSSCompute( task,
                             VSL_SS_COV|VSL_SS_COR,
                             VSL_SS_METHOD_FAST );
    CheckVslError(errcode);

    /***** Testing stat characteristics of mean and covariance matrix *****/
    dComputePvalsMean( dim, n, mean, a, (double*)C, pval_mean );
    dComputePvalsVariance( dim, n, (double*)cov, cov_storage, (double*)C,
                               (double*)pval_cov );
    dComputePvalsCovariance( dim, n, (double*)cov, cov_storage, (double*)C,
                             (double*)pval_cov );

    for(i = 0; i < dim; i++)
    {
        if (pval_mean[i] < P_THRESHOLD) errnums++;
        if (pval_cov[i][i] < P_THRESHOLD) errnums++;

        for(j = 0; j < i; j++)
        {
            if (pval_cov[i][j] < P_THRESHOLD) errnums++;
        }
    }

    /***** Printing results *****/
    printf("Task dimension :         %d\n", dim);
    printf("Number of observations : %d\n\n", n);

    /***** Print the exact mean, covariance and correlation matrices *****/
    printf("Exact means\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", a[i]);
    }

    printf("\n\nExact covariance matrix             ");
    printf("Exact correlation matrix\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+10lf ", C[i][j]);
        }

        printf("   ");

        for(j = 0; j < i; j++)
        {
            aTmp = C[i][i] * C[j][j];
            vdSqrt( 1, &aTmp, &rTmp );
            printf( "%+10lf ", C[i][j] / rTmp );
        }

        printf("%+10lf ", C[i][j]);

        for(j = i+1; j < dim; j++)
        {
            aTmp = C[i][i] * C[j][j];
            vdSqrt( 1, &aTmp, &rTmp );
            printf( "%+10lf ", C[i][j] / rTmp );
        }
        printf("\n");
    }

    /***** Print the computed mean, covariance and correlation matrices *****/
    printf("\nComputed means\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", mean[i]);
    }

    printf("\n\nComputed covariance matrix          ");
    printf("Computed correlation matrix\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+10lf ", cov[i][j]);
        }

        printf("   ");

        for(j = 0; j < dim; j++)
        {
            printf("%+10lf ", cor[i][j]);
        }
        printf("\n");
    }

    printf("\nP-values of the computed means\n");

    for(i = 0; i < dim; i++)
    {
        printf("%9lf ", pval_mean[i]);
    }

    printf("\n\nP-values of the computed covariance matrix\n");

    for(i = 0; i < dim; i++)
    {
        for(j = 0; j <= i; j++)
        {
            printf("%9lf ", pval_cov[i][j]);
        }
        printf("\n");
    }

    /***** Printing summary of the test *****/
    if (errnums == 0) {
        printf("\n\nMean and covariance estimates agree with theory\n");
    }
    else {
        printf("\n\nError: Mean and/or covariance estimates");
        printf(" disagree with theory\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
