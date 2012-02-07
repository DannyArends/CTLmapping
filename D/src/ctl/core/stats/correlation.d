/**********************************************************************
 * src/ctl/core/stats/correlation.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.core.stats.correlation;

import std.math;

//D routines for correlation analysis almost as fast as the RBlas version 
//But this allows us to write out to disk, using almost no additional RAM
double correlation(T)(T[] x, T[] y){ 
  double XiYi = 0;
  double Xi = 0;
  double Yi = 0;
  double XiP2 = 0;
  double YiP2 = 0;
  for(uint i = 0; i < x.length; i++){
    XiYi += x[i] * y[i];
    Xi += x[i]; 
    Yi += y[i];
    XiP2 += x[i] * x[i]; 
    YiP2 += y[i] * y[i];
  }
  double onedivn = 1.0/x.length;
  return (XiYi - (onedivn*Xi*Yi)) / (sqrt(XiP2 - onedivn * pow(Xi, 2)) * sqrt(YiP2 - onedivn * pow(Yi, 2)));
}
