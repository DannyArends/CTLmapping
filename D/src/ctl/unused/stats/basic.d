/******************************************************************//**
 * \file ctl/core/stats/basic.d
 * \brief Basic statistics routines
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.stats.basic;
 
import std.math, std.conv, std.string;

pure T doMean(T)(in T[] d){
  T mean = 0;
  for(size_t x = 0; x < d.length; x++){
    mean += (d[x] - mean) / (x + 1);
  }
  return mean;
}

T[] doMatrixMax(T)(in T[][] d){
  T[] m;
  m.length = d.length;
  for(size_t x=0; x < d.length; x++){ m[x] = doMax!T(d[x]); }
  return m;
}

pure T doMax(T)(in T[] d){
  T m = T.min;
  for(size_t x=0; x < d.length; x++){ if(m < d[x]){ m = d[x]; } }
  return m;
}

pure T doSum(T)(in T[] d){
  T s = 0;
  for(size_t x=0; x < d.length; x++){ s += d[x]; }
  return s;
}

pure real doSumOfSquares(T)(in T[] data){
  T mean = doMean(data);
  real sumofsquares = 0;
  for(size_t i = 0; i < data.length; i++){
    sumofsquares += pow((data[i]-mean),2);
  }
  return sumofsquares;
}

pure real doVariance(T)(T[] data){ return (doSumOfSquares!T(data)/(data.length-1)); }

pure real doVariance(T)(real sumofsquares,uint n){ return (sumofsquares/(n-1)); }

pure real doStandardDeviation(T)(T[] data){ return sqrt(doVariance!T(data)); }

unittest{
  import ctl.io.terminal;
  MSG("Unit test: %s", __FILE__);
  string test_fun;
  try{
    test_fun = "pure T doMean(T)(in T[] data)";
    assert(doMean([3.0, 4.0]) == 3.5,           "\n"~test_fun~" Test 1");
    assert(doMean([1.0, 2.0, 3.0, 4.0]) == 2.5, "\n"~test_fun~" Test 2");
    assert(doMean([1.0]) == 1.0,                "\n"~test_fun~" Test 3");
    assert(doMean([3, 4]) == 3,                 "\n"~test_fun~" Test 4");
    OK("Tests: %s",test_fun);

    test_fun = "pure T doMax(T)(in T[] d)";
    assert(doMax([3.0, 4.0]) == 4.0,           "\n"~test_fun~" Test 1");
    assert(doMax([1.0, 5.0, 3.0, 2.0]) == 5.0, "\n"~test_fun~" Test 2");
    assert(doMax([1.0]) == 1.0,                "\n"~test_fun~" Test 3");
    assert(doMax([3, 4]) == 4,                 "\n"~test_fun~" Test 4");
    OK("Tests: %s",test_fun);

    test_fun = "pure T doSum(T)(in T[] d)";
    assert(doSum([3.0, 4.0]) == 7.0,           "\n"~test_fun~" Test 1");
    assert(doSum([1.0, 2.3, 3.0, 1.0]) == 7.3, "\n"~test_fun~" Test 2");
    assert(doSum([1.0]) == 1.0,                "\n"~test_fun~" Test 3");
    assert(doSum([3, 4]) == 7,                 "\n"~test_fun~" Test 4");
    OK("Tests: %s",test_fun);

    MSG("Tested: %s",__FILE__);  
  }catch(Throwable e){
    string err = to!string(e).split("\n")[1];
    ERR("Reason: %s failed", err);
  }
}
