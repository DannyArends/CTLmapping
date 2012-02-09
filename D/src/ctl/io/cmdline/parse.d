/**********************************************************************
 * src/ctl/io/cmdline/parse.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.io.cmdline.parse;

import std.stdio;
import std.string;
import std.conv;

struct S(T){
  string name;
  string description;
  T      value;
}

struct Scontainer{
  S!bool[]   booleans;
  S!uint[]   integers;
  S!string[] strings;
}

struct CTLsettings{
  void load(){
    opt.booleans ~= S!bool("--help"         ,"Show the help file", false);
    opt.booleans ~= S!bool("--verbose"      ,"Verbose mode", true);
    opt.booleans ~= S!bool("--overwrite"    ,"Overwrite previous output files", false);
    opt.integers ~= S!uint("--nperms"       ,"Number of permutations", 100);
    opt.strings  ~= S!string("--phenotypes" ,"File containing phenotypes", "test/data/phenotypes.csv");
    opt.strings  ~= S!string("--genotypes"  ,"File containing genotypes", "test/data/genotypes.csv");
    opt.strings  ~= S!string("--format"     ,"File format", "csv");
  }
  
  bool displayHelp(){
    if(getBool("--help")){
      foreach(S!bool s; opt.booleans){ writeln("  " ~ s.name[0] ~ "[" ~ s.name[1..3] ~"]" ~ s.name[3..$] ~ " - " ~ s.description); }
      foreach(S!uint s; opt.integers){ writeln("  " ~ s.name[0] ~ "[" ~ s.name[1..3] ~"]" ~ s.name[3..$] ~ " - " ~ s.description); }
      foreach(S!string s; opt.strings){ writeln("  " ~ s.name[0] ~ "[" ~ s.name[1..3] ~"]" ~ s.name[3..$] ~ " - " ~ s.description); }
      return true;
    }
    return false;
  }
  
  uint getInt(string name){
    foreach(int cnt, S!uint s; opt.integers){
      if(name==s.name){
        return s.value;
      }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return -1;
  }
  
  bool getBool(string name){
    foreach(int cnt, S!bool s; opt.booleans){
      if(name==s.name){
        return s.value;
      }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return false;
  }
  
  string getString(string name){
    foreach(int cnt, S!string s; opt.strings){
      if(name==s.name){
        return s.value;
      }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return "";
  }
  Scontainer opt;
}

CTLsettings parseCmd(string[] args){
  CTLsettings settings;
  settings.load();
  if(args.length > 1){
    foreach(string arg; args[1..$]){
      bool known = false;
      auto isloc = arg.indexOf("=");
      if(isloc > 1){
        if(arg.length >= 4){  // -p=a -g=t.txt <- minimum of 4
          foreach(int cnt, S!string s; settings.opt.strings){
            if(arg[1]==s.name[2]){
              try{
              settings.opt.strings[cnt].value = arg[(isloc+1)..$];
              known = true;
              writeln(" - CMD line: ", s.description, " set to ", settings.opt.strings[cnt].value);
              }catch(Throwable e){ }
            }
          }
          foreach(int cnt, S!uint s; settings.opt.integers){
            if(arg[1]==s.name[2]){
              try{
              settings.opt.integers[cnt].value = to!int(arg[(isloc+1)..$]);
              known = true;
              writeln(" - CMD line: ", s.description, " set to ", settings.opt.integers[cnt].value);
              }catch(Throwable e){ }
            }
          }
        }
      }else{ // Boolean perhaps ?
        if(arg.length == 2){ // -v -h <- exactly two long
          foreach(int cnt, S!bool s; settings.opt.booleans){
            if(arg[1]==s.name[2]){
              known = true;
              writeln(" - CMD line: " ~ s.description);
              settings.opt.booleans[cnt].value = !settings.opt.booleans[cnt].value;
            }
          }
        }
      }
      if(!known) writeln("Unknown argument: ",arg);
    }
  }
  return settings;  
}
