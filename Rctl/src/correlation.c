/******************************************************************//**
 * \file Rctl/src/correlation.c
 * \brief Implementation of functions related to correlation
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "correlation.h"

/* Calculates pearsons correlation between x and y (Ranked input for non-parametric testing) */
double correlation(const double* x, const double* y, size_t dim, bool verbose){
  size_t i;
  double onedivn = 1.0 / dim;
  KahanAccumulator XiYi = createAccumulator();
  KahanAccumulator Xi   = createAccumulator();
  KahanAccumulator Yi   = createAccumulator();
  KahanAccumulator XiP2 = createAccumulator();
  KahanAccumulator YiP2 = createAccumulator();

  for(i = 0; i < dim; i++){
    if(x[i] != MISSING && y[i] != MISSING){
      XiYi = KahanSum(XiYi, x[i] * y[i]);
      Xi   = KahanSum(Xi, x[i]);
      Yi   = KahanSum(Yi, y[i]);
      XiP2 = KahanSum(XiP2, x[i] * x[i]);
      YiP2 = KahanSum(YiP2, y[i] * y[i]);
    }
  }
  double nom = (XiYi.sum - (onedivn*Xi.sum*Yi.sum));
  double denom = sqrt(XiP2.sum - onedivn * pow(Xi.sum, 2.0)) * sqrt(YiP2.sum - onedivn * pow(Yi.sum, 2.0));
  double cor = nom/denom;
  if(isNaN(cor) || isinf(cor) || cor < -1 || cor > 1){ 
    err("Correlation '%.4f' not in range [-1, 1]\n", cor); 
  }
  return(cor);
}

/* Calculate the chi square test statistic based on N seggregating correlations */
double chiSQ(size_t nr, double* r, int* nsamples){
  size_t i;
  size_t denom = 0;
  KahanAccumulator sumOfSquares = createAccumulator();
  KahanAccumulator squaresOfSum = createAccumulator();

  for(i = 0; i < nr; i++){
    sumOfSquares = KahanSum(sumOfSquares, (nsamples[i]-3) * pow(zscore(r[i]), 2.0));
    squaresOfSum = KahanSum(squaresOfSum, (nsamples[i]-3) * zscore(r[i]));
    denom  += (nsamples[i]-3);
  }
  if(denom == 0) info("Divide by 0 groups too small");
  return(sumOfSquares.sum - (pow(squaresOfSum.sum, 2.0) / denom));
}

/* Transforms a chi square critical value (Cv) to a p-value */
double chiSQtoP(int Dof, double Cv){
  if(Cv <= 0 || Dof < 1) return 1.0;
  return dchisq(Cv,(double)Dof, 0);
}

