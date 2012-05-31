/******************************************************************//**
 * \file ctl/core/stats/correlation.d
 * \brief Basic correlation implementation
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.stats.correlation;

import std.math;

//D routines for correlation analysis almost as fast as the RBlas version 
//But this allows us to write out to disk, using almost no additional RAM
double correlation(T)(in T[] x, in T[] y){
  assert(x.length == y.length);
  double XiYi = 0;
  double Xi = 0;
  double Yi = 0;
  double XiP2 = 0;
  double YiP2 = 0;
  for(size_t i = 0; i < x.length; i++){
    XiYi += x[i] * y[i];
    Xi += x[i]; 
    Yi += y[i];
    XiP2 += x[i] * x[i]; 
    YiP2 += y[i] * y[i];
  }
  double onedivn = 1.0/x.length;
  return (XiYi - (onedivn*Xi*Yi)) / (sqrt(XiP2 - onedivn * pow(Xi, 2)) * sqrt(YiP2 - onedivn * pow(Yi, 2)));
}
