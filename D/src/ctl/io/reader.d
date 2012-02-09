/**********************************************************************
 * src/ctl/io/reader.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
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
  abstract int[][] loadgenotypes(string filename);
}
