/******************************************************************//**
 * \file ctl/core/array/matrix.d
 * \brief Matrix and vector functions
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.array.matrix;

import std.stdio, std.conv, std.math, core.memory;

T[][] absmatrix(T)(in T[][] i){
  T[][] m = newmatrix!T(i.length, i[0].length);
  for(size_t r=0;r<i.length;r++){
    for(size_t c=0;c<i[0].length;c++){
      m[r][c] = abs(i[r][c]);
    }
  }
  return m;
}

T[] unlist(T)(in T[][] i){
  T[] m;
  for(size_t r=0;r<i.length;r++){ m ~= i[r]; }
  return m;
}

T[][] translate(T)(in T[][] i){
  T[][] m = newmatrix!T(i[0].length,i.length);
  for(uint r=0;r<i.length;r++){
    for(uint c=0;c<i[0].length;c++){
      m[c][r] = i[r][c];
    }
  }
  return m;
}

T[][] newmatrix(T)(size_t rows, size_t cols, T value = T.init){
  T[][] m;
  m.length = rows;
  for(size_t i = 0; i < rows; i++){
    m[i] = newvector!T(cols,value);
  }
  return m;
}

void printmatrix(T)(T[][] m, string sep = " ") {
  for(int r=0; r < m.length; r++){
    for(int c=0; c < m[r].length; c++){
      write(to!string(m[r][c]), sep);
    }
    writeln();
  }
}

T[] newvector(T)(size_t dim, T value = T.init) {
  T[] v;
  v.length = dim;
  for(int e=0; e<dim; e++){ v[e] = value; }
  return v;
}

void freevector(T)(ref T[] v) {
  GC.removeRange(cast(void*)v);
  GC.free(cast(void*)v);
}

void freematrix(T)(T[][] m) {
  for(size_t i=0; i < m.length; i++) {
    if(m[i].length > 0) freevector!T(m[i]);
  }
  GC.removeRange(cast(void*)m);
  GC.free(cast(void*)m);
}

T[] stringvectortotype(T)(string[] entities){
  T[] rowleveldata;
  rowleveldata.length = entities.length;
  for(size_t cnt = 0;cnt < entities.length; cnt++){
    try{
      rowleveldata[cnt] = to!T(entities[cnt]);
    }catch(Throwable e){
      rowleveldata[cnt] = T.max;
    }
  }
  return rowleveldata;
}
