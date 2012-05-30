/******************************************************************//**
 * \file ctl/io/terminal.d
 * \brief Functions to modify the terminal
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written feb, 2011<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.io.terminal;
private import std.stdio;

static{
  CColor fg = CColor.White;
  CColor bg = CColor.Black;
}

version(Windows){ 
  import std.c.windows.windows;

  enum CColor : ushort{ 
    Black  = 0 , DarkBlue   = 1 , DarkGreen = 2 , DarkAzure = 3 , DarkRed = 4,
    Purple = 5 , DarkYellow = 6 , Silver    = 7 , Gray      = 8 , Blue    = 9,
    Green  = 10, Aqua       = 11, Red       = 12, Magenta   = 13, Yellow  = 14,
    White  = 15, Default    = 256}
  
  static{ extern(C) HANDLE hConsole = null; }

  static this(){
    hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO info;
    GetConsoleScreenBufferInfo( hConsole, &info );
    bg = cast(CColor)(info.wAttributes & (0b11110000));
    fg = cast(CColor)(info.wAttributes & (0b00001111));
  }

  private ushort asCColor(CColor fg, CColor bg){ return cast(ushort)(fg | bg << 4); }

  private void setConsoleForeground(CColor color = CColor.White){
    if(color != fg){ SetConsoleTextAttribute(hConsole, asCColor(color, bg)); }
    fg = color;
  }

  private void setConsoleBackground(CColor color = CColor.Black){   
    if(color != bg){ SetConsoleTextAttribute(hConsole, asCColor(fg, color)); }
    bg = color;
  }  

}else version(Posix){
    
  enum CColor : ushort{
    Black = 30, Red  = 31, Green = 32, Yellow  = 33, Blue = 34,
    Pink  = 35, Aqua = 36, White = 37, Default = 0 }
    
  private void setConsoleForeground(CColor color){
    if(color == CColor.Default){
      writef("\033[0m");
      fg = CColor.Default;
      setConsoleBackgroundColor(bg);
    }else{
      writef("\033[0;%dm", cast(int)color);
      fg = color;
    }
  }

  void setConsoleBackground(CColor color){
    if(color == CColor.Default){
      writef("\033[0m");
      bg = CColor.Default;
      setConsoleFontColor(fg);
    }else{
      writef("\033[0;%dm", cast(int)color + 10);
      bg = color;
    }
  }  
}

void setCColor(CColor fore = CColor.White, CColor back = CColor.White.Black){
  setConsoleForeground(fore);
  setConsoleBackground(back);
}

void addWhite(string name){
  assert(name.length < 6, "The maximum name length of GenOutput is 6");
  for(size_t x=0; x < (6-name.length);x++){ write(" "); }
}

template GenOutput(string name, string color, string bcolor = "Black"){
  const char[] GenOutput = "void w"~ name ~ "(A...)(A a){" ~
  "setCColor(CColor."~color~",CColor."~bcolor~"); writef(\"["~name~"]\"); std.stdio.stdout.flush();" ~
  "setCColor(); addWhite(\""~name~"\"); std.stdio.stdout.flush();" ~
  "setCColor(); writefln(a);}";
}

mixin(GenOutput!("CTL", "Green", "Black"));
alias wCTL MSG;
mixin(GenOutput!("ERROR", "Red", "Black"));
alias wERROR ERR;
mixin(GenOutput!("WARN", "Yellow", "Black"));
alias wWARN WARN;
mixin(GenOutput!("DEBUG", "White", "Black"));
alias wDEBUG DBG;

CColor getConsoleForeground(){ return fg; }
CColor getConsoleBackground(){ return bg; }
