/******************************************************************//**
 * \file ctl/io/csv/write.d
 * \brief Write output in CSV format
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.io.csv.write;

import std.stdio, std.string, std.file, std.conv;
import ctl.core.array.matrix, ctl.io.terminal;

void writeFile(T)(T[][] m, string filename, string[] rownames, bool overwrite = false, bool verbose = false){
  if(exists(filename)){
    if(overwrite){
      remove(filename);
    }else{ return; }
  }
  try{
    auto fp = new File(filename,"wb");
    string      buffer;
    for (int r=0; r<m.length; r++) {
      if(rownames != null)  fp.write(rownames[r],"\t");
      for (int c=0; c<m[r].length; c++) {
        if(c!=0) fp.write("\t");
        fp.write(to!string(m[r][c]));
      }
      fp.writeln();
    }
    fp.close();
  }catch(Throwable e){
    writefln("File %s write exception: %s", filename,e);
  }
}


void addToFile(T)(T[][] m, string filename){
  try{
    MSG("addToFile");
    auto fp = new File(filename,"a");
    string      buffer;
    for (int r=0; r<m.length; r++) {
      for (int c=0; c<m[r].length; c++) {
        if(c!=0) fp.write("\t");
        fp.write(to!string(m[r][c]));
      }
      fp.writeln();
    }
    fp.close();
  }catch(Throwable e){
    writefln("File %s write exception: %s", filename,e);
  }
}
