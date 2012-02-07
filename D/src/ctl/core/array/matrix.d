/**********************************************************************
 * src/ctl/core/array/matrix.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Feb, 2012
 * first written May, 2011
 **********************************************************************/
module ctl.core.array.matrix;

import std.stdio;
import std.conv;
import std.math;

T[][] absmatrix(T)(T[][] i){
  T[][] m = newmatrix!T(i.length,i[0].length);
  for(uint r=0;r<i.length;r++){
    for(uint c=0;c<i[0].length;c++){
      m[r][c] = abs(i[r][c]);
    }
  }
  return m;
}

T[] unlist(T)(T[][] i){
  T[] v = newvector!T(i.length*i[0].length);
  for(uint r=0;r<i.length;r++){
    v[r*i.length..(r*i.length)+i[0].length] = i[r];
  }
  return v;
}

T[][] translate(T)(T[][] i){
  T[][] m = newmatrix!T(i[0].length,i.length);
  for(uint r=0;r<i.length;r++){
    for(uint c=0;c<i[0].length;c++){
      m[c][r] = i[r][c];
    }
  }
  return m;
}

T[][] newmatrix(T)(size_t rows, size_t cols) {
  T[][] m;
  m.length=rows;
  if(m is null){
    writeln("Not enough memory for new matrix");
  }
  for(size_t i=0; i<rows; i++) {
    m[i].length= cols;
  }
  return m;
}

void printmatrix(T)(T[][] m) {
  for (int r=0; r<m.length; r++) {
    for (int c=0; c<m[r].length; c++) {
      write(to!string(m[r][c])," ");
    }
    writeln();
  }
}

T[] newvector(T)(size_t dim) {
  T[] v;
  v.length = dim;
  if(v is null){
    writeln("Not enough memory for new vector of dimension %d",(dim+1));
  }
  return v;
}

T[] stringvectortotype(T)(string[] entities){
  T[] rowleveldata;
  for(auto e=0;e < entities.length; e++){
    try{
      rowleveldata ~= to!T(entities[e]);
    }catch(Throwable e){
      rowleveldata ~= to!T(0);
    }
  }
  return rowleveldata;
}
