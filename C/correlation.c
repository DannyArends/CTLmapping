#include "correlation.h"
 
double correlation(double* x, double* y, size_t dim){
  double XiYi = 0.0;
  double Xi   = 0.0;
  double Yi   = 0.0;
  double XiP2 = 0.0;
  double YiP2 = 0.0;
  size_t i    = 0;
  for(i = 0; i < dim; i++){
    XiYi += x[i] * y[i];
    Xi   += x[i];
    Yi   += y[i];
    XiP2 += x[i] * x[i]; 
    YiP2 += y[i] * y[i];
  }
  double onedivn = 1.0 / dim;
  return (XiYi - (onedivn*Xi*Yi)) / (sqrt(XiP2 - onedivn * pow(Xi, 2)) * sqrt(YiP2 - onedivn * pow(Yi, 2)));
}
