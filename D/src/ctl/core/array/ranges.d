/**********************************************************************
 * src/ctl/core/array/ranges.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Feb, 2012
 * first written May, 2011
 **********************************************************************/
module ctl.core.array.ranges;

import std.math;

pure uint[] dorange(int start, size_t length){
  uint array[];
  array.length = length-1;
  for(uint i = 0; i < (length-1); i++){
    array ~= start+i;
  }
  return array;
}

pure T[] doarray(T)(int length,T value){
  T array[];
  array.reserve(length-1);
  for(int i = 0; i < (length-1); i++){
   array ~= value;
  }
  return array;
}
