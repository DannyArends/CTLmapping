/******************************************************************//**
 * \file ctl/core/stats/correlation.d
 * \brief Basic correlation implementation
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified Sep, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.stats.correlation;

import std.math, std.stdio, std.traits, std.string;

//D routines for correlation analysis almost as fast as the RBlas version 
//But this allows us to write out to disk, using almost no additional RAM
T correlation(T = double)(in T[] x, in T[] y)
in{
  assert(isFloatingPoint!T, "T must be a FloatingPoint type");
  assert(x.length == y.length, "x.length must match y.length");
}
out(result){
  assert((round(result) >= -1) && (round(result) <= 1), format("Correlation out of bounds: -1 <= result <= 1 = %s",result));
}
body{
  T XiYi = 0;
  T Xi = 0;
  T Yi = 0;
  T XiP2 = 0;
  T YiP2 = 0;
  for(size_t i = 0; i < x.length; i++){
    XiYi += x[i] * y[i];
    Xi   += x[i]; 
    Yi   += y[i];
    XiP2 += x[i] * x[i]; 
    YiP2 += y[i] * y[i];
  }
  T onedivn = 1.0/x.length;
  return (XiYi - (onedivn*Xi*Yi)) / (sqrt(XiP2 - onedivn * pow(Xi, 2)) * sqrt(YiP2 - onedivn * pow(Yi, 2)));
}
