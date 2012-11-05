#include "correlation.h"
 
double correlation(double* x, double* y, size_t dim){
  double XiYi = 0.0;
  double Xi   = 0.0;
  double Yi   = 0.0;
  double XiP2 = 0.0;
  double YiP2 = 0.0;
  size_t i    = 0;
  for(i = 0; i < dim; i++){
    if(x[i] != -999 && y[i] != -999){
      XiYi += x[i] * y[i];
      Xi   += x[i];
      Yi   += y[i];
      XiP2 += x[i] * x[i];
      YiP2 += y[i] * y[i];
    }
  }
  double onedivn = 1.0 / dim;
  return (XiYi - (onedivn*Xi*Yi)) / (sqrt(XiP2 - onedivn * pow(Xi, 2)) * sqrt(YiP2 - onedivn * pow(Yi, 2)));
}

int isNaN(double d){ 
  return(d != d); 
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
