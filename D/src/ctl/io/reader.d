/******************************************************************//**
 * \file ctl/io/reader.d
 * \brief Reader interface to support multiple file formats
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.io.reader;

import ctl.io.csv.parse;
version(QTAB){
  import ctl.io.qtab.wrapper;
}
import ctl.io.cmdline.parse;

Reader initialize(CTLsettings settings){
  switch(settings.getString("--format")){
version(QTAB){
    case "qtab":
     return new QTABreader();
    break;
}
    default:
    break;
  }
  return new CSVreader();
}


abstract class Reader{
  abstract double[][] loadphenotypes(string filename);
  abstract string[] loadphenonames(string filename);
  abstract int[][] loadgenotypes(string filename);
}
