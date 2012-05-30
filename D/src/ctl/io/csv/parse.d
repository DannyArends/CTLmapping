/******************************************************************//**
 * \file ctl/io/csv/parse.d
 * \brief Parse input in CSV format
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.io.csv.parse;

import std.stdio, std.string, std.file, std.conv;
import ctl.io.reader, ctl.io.terminal, ctl.core.array.matrix;

class CSVreader : Reader{
  double[][] loadphenotypes(string filename = "phenotypes.csv"){
    return parseFile!double(filename);
  }

  int[][] loadgenotypes(string filename = "genotypes.csv"){
    return parseFile!int(filename);
  }
}

struct D{
  string  name;
  char    chr = '0';
  double  loc = 0.0;
  
  this(string[] d){
    name = d[0];
    try{
    if(d[1] != "") chr  = to!char(d[1]);
    if(d[2] != "") loc  = to!double(d[2]);
    }catch(Throwable t){ }
  }
}

T[][] parseFile(T)(string filename, bool verbose = false ,bool hasRowHeader= true){
  T[][] data;
  D[]   descr;
  if(!exists(filename) || !isFile(filename)){
    ERR("No such file %s",filename);
  }else{
    try{
      auto fp = new File(filename,"rb");
      string      buffer;
      while(fp.readln(buffer)){
        buffer = chomp(buffer);
        string[] splitted = buffer.split("\t");
        if(hasRowHeader){
          if(splitted.length > 3){
            descr ~= D(splitted[0..2]);
            data ~= stringvectortotype!T(splitted[3..$]);
          }
        }else{
          data ~= stringvectortotype!T(splitted[0..$]);        
        }
      }
      fp.close();
      if(verbose) MSG("Parsed %s imports from file: %s",data.length, filename);
    }catch(Throwable e){ ERR("File %s read exception", filename); }
  }
  return data;
}
