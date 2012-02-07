/**********************************************************************
 * src/ctl/io/csv/write.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.io.csv.write;

import std.stdio;
import std.string;
import std.file;
import std.conv;

import ctl.core.array.matrix;

void writeFile(T)(T[][] m, string filename, bool verbose = false){
  if(exists(filename)){
    writefln("Already exists %s",filename);
  }else{
    try{
      auto fp = new File(filename,"wb");
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
      writefln("File %s read exception: %s", filename,e);
    }
  }
}
