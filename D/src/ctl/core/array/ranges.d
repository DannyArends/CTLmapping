/******************************************************************//**
 * \file ctl/core/array/ranges.d
 * \brief Range functions
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.array.ranges;

import std.math;

pure T max(T)(T v, T v2){if(v < v2){return v2;}else{return v;}}

pure T max(T)(in T[] r){
  assert(r.length >= 0);
  T best = T.min;
  foreach(e; r){
    if(e > best) best = e;
  }
  return best;
}

pure size_t[] dorange(int start, size_t length){
  size_t array[];
  for(size_t i = 0; i < length; i++){ array ~= start+i; }
  return array;
}
