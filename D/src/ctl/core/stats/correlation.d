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

import std.math, std.stdio, std.math, std.mathspecial, std.conv, std.string;

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

double studentT(T)(in T[] x, in T[] y){
  double corr = correlation(x,y);
  double df = x.length-2.0;
  return studentTC(corr,df);
}

double studentTC(double corr, double df){
  double t = abs(corr * sqrt(df / (1 - corr*corr)));
  writefln("Cor:%s", corr);
  writefln("T:%s",  t);
  return t;
}

double pvalue(T)(in T[] x, in T[] y){
  double t = studentT(x,y);
  double df = x.length-2.0;
  double k = gamma((df+1) / 2.0);
  double k2 = k / (sqrt(PI*df) * gamma(df / 2.0));
  double p = pow(k2 * (pow(1+t,2)/df), -((df+1)/2.0));
  writefln("DF:%s\nP:%s", df, p);
  return p;
}

double zscore(T)(in T[] x, in T[] y){
  double p = pvalue(x,y);
  writefln("Z:%s",normalDistributionInverse(p));
  return normalDistributionInverse(p);
}

unittest{
  import ctl.io.terminal;
  MSG("Unit test: %s", __FILE__);
  string test_fun;
  try{
    test_fun = "double correlation(T)(in T[] x, in T[] y)";
    for(double x = -1;x <= 1;x+=0.1){
      writefln("%s T: %s",x,studentTC(x,10));    
    }
    writefln("Z: %s",zscore([1,3,4,7,8,10,15],[1,3,4,7,8,10,15]));
    writefln("Z: %s",zscore([1,3,4,7,8,10,15,20],[10,3,5,7,3,1,12,4]));
    assert(correlation([1,3,4,7,8,10,15],[1,3,4,7,8,10,15]),    "\n"~test_fun~" Test 1");
    OK("Tests: %s",test_fun);

    MSG("Tested: %s",__FILE__);  
  }catch(Throwable e){
    string err = to!string(e).split("\n")[1];
    ERR("Reason: %s failed", err);
  }
}
