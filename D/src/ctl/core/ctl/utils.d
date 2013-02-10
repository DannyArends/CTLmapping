/******************************************************************//**
 * \file ctl/core/ctl/utils.d
 * \brief Utility classes for CTL mapping
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.ctl.utils;

import std.stdio, std.conv, std.string, std.c.stdlib;
import ctl.core.array.matrix;

/* Write an error string to stderr */
void error(in string s){
  stderr.writeln();
  stderr.writefln("-Error: %s\n", s);
}

/* Abort with error code, default: -1 */
void abort(in string s, int exitcode = -1){
  error(s);
  exit(exitcode);
}

