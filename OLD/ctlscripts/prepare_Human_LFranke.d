import std.stdio,std.conv, std.file, std.string, std.array, std.math;


T[] stringstotype(T)(in string[] entities){
  T[] rowleveldata;
  for(size_t e=0; e < entities.length; e++){
    try{
      rowleveldata ~= to!T(entities[e]);
    }catch(Throwable e){
      rowleveldata ~= T.init;
    }
  }
  return rowleveldata;
}

string[] myround(double[] i, int d = 3){
  string[] rounds;
  for(size_t e=0; e < i.length; e++){
    rounds ~= to!string(round(i[e]* pow(10,d)) / to!double(pow(10,d)));
  }
  return rounds;
}

void main(string[] args){
  string fn = "input/BloodSATVATLiverMuscleHT12ProbesCentered.txt.TriTyperFormat.txt";
  if(exists(fn) && isFile(fn)){
    string buffer;
    int cnt = 0;
    string fnout = "out_RAW.txt";
    auto f = new File(fn,"rb");
    auto fout = new File(fnout,"w");
    f.readln(buffer);
    string[] header = chomp(buffer).split("\t");
    fout.writeln(header[0],"\t", join(header[9..$],"\t"));
    while(f.readln(buffer)){
      string[] items = chomp(buffer).split("\t");
      string[] values = myround(stringstotype!double(items[9..$]));
      fout.writeln(items[0],"\t", join(values, "\t"));
      cnt++;
      if(cnt % 1000 == 0) writefln("Processed %s lines",cnt);
    }
    writefln("Done after: %s lines",cnt);
  }
}
