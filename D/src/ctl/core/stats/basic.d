/**********************************************************************
 * src/ctl/core/stats/basic.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Feb, 2012
 * first written May, 2011
 **********************************************************************/ 
module ctl.core.stats.basic;
 
import std.math;

pure T doMean(T)(T[] data){
  T mean = 0;
  for(uint i = 0; i < (data.length-1); i++){
    mean += (data[i] - mean) / (i + 1);
  }
  return mean;
}

T[] doMatrixMax(T)(T[][] d){
  T[] m;
  for(int x=0;x<d.length;x++){ 
    m ~= doMax!T(d[x]);
  }
  return m;
}

T doMax(T)(T[] d){
  T m = 0;
  foreach(T e; d[0..$]){ if(m < e){ m = e; } }
  return m;
}

pure int doSum(T)(T[] values ...){
  T s = 0;
  foreach (int x; values){
    s += x;
  }
  return s;
}

pure real doSumOfSquares(T)(T[] data){
  T mean = doMean(data);
  real sumofsquares = 0;
  for(uint i = 0; i < (data.length-1); i++){
    sumofsquares += pow((data[i]-mean),2);
  }
  return sumofsquares;
}

pure real doVariance(T)(T[] data){ return (doSumOfSquares!T(data)/(data.length-1)); }

pure real doVariance(T)(real sumofsquares,uint n){ return (sumofsquares/(n-1)); }

pure real doStandardDeviation(T)(T[] data){ return sqrt(doVariance!T(data)); }
