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

import std.stdio, std.string, std.file, std.conv, std.datetime;
import ctl.io.reader, ctl.core.array.matrix, ctl.core.ctl.utils;

class CSVreader : Reader{
  override double[][] loadphenotypes(string filename = "phenotypes.csv"){
    return parseFile!double(filename,false,true,0);
  }
  
  override string[] loadphenonames(string filename = "phenotypes.csv"){
    return parseNames(filename);
  }
  
  override int[][] loadgenotypes(string filename = "genotypes.csv"){
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

string[] parseNames(string filename){
  string[] data;
  if(!exists(filename) || !isFile(filename)){
    abort(xformat("No such file %s", filename));
  }else{
    try{
      SysTime stime = Clock.currTime();
      string[] content = readText(filename).split("\n");
      foreach(string buffer; content){
        string[] splitted = chomp(buffer).split("\t");
        if(splitted.length > 0){ data ~= splitted[0]; }
      }
    }catch(Throwable e){ abort(xformat("File %s read exception %s", filename, e)); }
  }
  return data;
}


T[][] parseFile(T)(string filename, bool verbose = false ,bool hasRowHeader = true, T nullval = -999){
  T[][] data;
  if(!exists(filename) || !isFile(filename)){
    abort(xformat("No such file %s", filename));
  }else{
    try{
      SysTime stime = Clock.currTime();
      auto f = new File(filename,"rb");
      uint filesize = cast(uint)getSize(filename);
      ubyte[] inputbuffer = new ubyte[filesize];
      f.rawRead(inputbuffer);
      f.close();
      string[] content = (cast(string)inputbuffer).split("\n");
      string buffer;
      for(size_t cnt=0; cnt < content.length;cnt++){
        buffer = content[cnt];
        string[] splitted = chomp(buffer).split("\t");
        if(splitted.length > 0){
        if(hasRowHeader){
          data ~= stringvectortotype!T(splitted[3..$],nullval);
        }else{
          data ~= stringvectortotype!T(splitted[0..$],nullval);
        }
        }
        freevector(splitted);
        //if((Clock.currTime()-stime).total!"seconds"() > 0 && (Clock.currTime()-stime).total!"seconds"() % 10 == 0) MSG("At: %s", data.length);
      }
      freevector(content);
      freevector(inputbuffer);
      if(verbose) writefln("Parsed %s imports from file: %s",data.length, filename);
      writefln("Loading took: (%s msecs)",(Clock.currTime()-stime).total!"msecs"());
    }catch(Throwable e){ abort(xformat("File %s read exception %s", filename, e)); }
  }
  return data;
}
