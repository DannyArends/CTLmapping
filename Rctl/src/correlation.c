#include "correlation.h"

/* Calculates pearsons correlation between x and y (Ranked input for non-parametric testing) */
double correlation(double* x, double* y, size_t dim){
  size_t i;
  double onedivn = 1.0 / dim;
  KahanAccumulator XiYi = createAccumulator();
  KahanAccumulator Xi   = createAccumulator();
  KahanAccumulator Yi   = createAccumulator();
  KahanAccumulator XiP2 = createAccumulator();
  KahanAccumulator YiP2 = createAccumulator();

  for(i = 0; i < dim; i++){
    if(x[i] != -999 && y[i] != -999){
      XiYi = KahanSum(XiYi, x[i] * y[i]);
      Xi   = KahanSum(Xi, x[i]);
      Yi   = KahanSum(Yi, y[i]);
      XiP2 = KahanSum(XiP2, x[i] * x[i]);
      YiP2 = KahanSum(YiP2, y[i] * y[i]);
    }
  }
  double nom = (XiYi.sum - (onedivn*Xi.sum*Yi.sum));
  double denom = sqrt(XiP2.sum - onedivn * pow(Xi.sum, 2.0)) * sqrt(YiP2.sum - onedivn * pow(Yi.sum, 2.0));
  return nom / denom;
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

/* Inverse Gamma function */
double igf(double S, double Z){
  if(Z < 0.0) return 0.0;

  long double Sc = (1.0 / S);
  Sc *= pow(Z, S);
  Sc *= exp(-Z);
 
  KahanAccumulator sum = createAccumulator();//= 1.0;
  sum.sum = 1.0;
  long double nom = 1.0;
  long double denom = 1.0;
  size_t i;
  for(i = 0; i < 1000; i++){
    nom *= Z;
	  S++;
	  denom *= S;
	  sum = KahanSum(sum, (nom / denom));
  }
  return sum.sum * Sc;
}

/* Transforms a chi square critical value (Cv) to a p-value */
double chiSQtoP(int Dof, double Cv){
  if(Cv <= 0 || Dof < 1) return 1.0;
  return dchisq(Cv,(double)Dof, 0);
}

