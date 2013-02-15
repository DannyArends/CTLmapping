#include "correlation.h"
 
double correlation(double* x, double* y, size_t dim){
  KahanAccumulator XiYi = createAccumulator();
  KahanAccumulator Xi   = createAccumulator();
  KahanAccumulator Yi   = createAccumulator();
  KahanAccumulator XiP2 = createAccumulator();
  KahanAccumulator YiP2 = createAccumulator();
  size_t i    = 0;
  for(i = 0; i < dim; i++){
    if(x[i] != -999 && y[i] != -999){
      XiYi = KahanSum(XiYi, x[i] * y[i]);
      Xi   = KahanSum(Xi, x[i]);
      Yi   = KahanSum(Yi, y[i]);
      XiP2 = KahanSum(XiP2, x[i] * x[i]);
      YiP2 = KahanSum(YiP2, y[i] * y[i]);
    }
  }
  double onedivn = 1.0 / dim;
  return (XiYi.sum - (onedivn*Xi.sum*Yi.sum)) / (sqrt(XiP2.sum - onedivn * pow(Xi.sum, 2.0)) * sqrt(YiP2.sum - onedivn * pow(Yi.sum, 2.0)));
}

double chiSQ(size_t nr, double* r, int* nsamples){
  size_t i;
  KahanAccumulator sumOfSquares = createAccumulator();
  KahanAccumulator squaresOfSum = createAccumulator();
  size_t denominator = 0;
  for(i = 0; i < nr; i++){
    sumOfSquares = KahanSum(sumOfSquares, (nsamples[i]-3) * pow(zscore(r[i]), 2.0));
    squaresOfSum = KahanSum(squaresOfSum, (nsamples[i]-3) * zscore(r[i]));
    denominator  += (nsamples[i]-3);
  }
  return(sumOfSquares.sum - (pow(squaresOfSum.sum, 2.0) / denominator));
}

double tstat(double cor, int dim){
  double n = (1.0 - pow(cor,2));
  double d = cor * sqrt((dim-2) / n);
  if(n!=0.0 && !isNaN(d)){
    return(d);
  }else{
    info("NaN produced from %f, %i\n", cor, dim);
    return(0.0);
  }
}

