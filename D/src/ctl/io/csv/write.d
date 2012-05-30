/**********************************************************************
 * src/ctl/io/csv/write.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.io.csv.write;

import std.stdio, std.string, std.file, std.conv;
import ctl.core.array.matrix, ctl.io.terminal;

void writeFile(T)(T[][] m, string filename, bool overwrite = false, bool verbose = false){
  if(exists(filename)){
    if(overwrite){
      remove(filename);
    }else{ return; }
  }
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
