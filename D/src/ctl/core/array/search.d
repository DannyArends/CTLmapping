/**********************************************************************
 * src/ctl/core/array/search.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Feb, 2012
 * first written May, 2011
 **********************************************************************/
module ctl.core.array.search;

import std.math, std.conv, std.string;

pure size_t getIndex(T)(T[] haystack, T needle){
  return searchArrayBinary!T(haystack, needle);
}

pure bool searchArray(T)(T[] haystack, T needle){
  foreach(T s; haystack){
    if(s==needle) return true;
  }
  return false;
}

pure size_t searchArrayBinary(T)(T[] haystack, T needle){
  size_t first = 0;
  size_t last = (haystack.length-1);
  while(first <= last){
    if(last==first) return first;
    size_t mid = (first + last) / 2;
    if(needle > haystack[mid]){
      first = mid + 1;
    }else if(needle < haystack[mid]){
      last = mid - 1;
    }else{
      return mid;
    }
  }
  return last;
}

unittest{
  import ctl.io.terminal;
  MSG("Unit test: %s", __FILE__);
  string test_fun;
  try{
    test_fun = "pure size_t getIndex(T)(T[] haystack, T needle)";
    assert(getIndex([3,5],4) == 1,                 "\n"~test_fun~" Test 1");
    assert(getIndex([1,3,4,7,8,10,15],8)  == 4,    "\n"~test_fun~" Test 2");
    assert(getIndex([1,2,4,7,8,10,15],3)  == 2,    "\n"~test_fun~" Test 3");
    assert(getIndex([1,2,4,7,8,10,15],9)  == 4,    "\n"~test_fun~" Test 4");
    assert(getIndex([1,2,4,7,8,10,15],12) == 6,    "\n"~test_fun~" Test 5");
    assert(getIndex([1,2,4,7,8,10,15],29) == 6,    "\n"~test_fun~" Test 6");
    OK("Tests: %s",test_fun);

    MSG("Tested: %s",__FILE__);  
  }catch(Throwable e){
    string err = to!string(e).split("\n")[1];
    ERR("Reason: %s failed", err);
  }
}

