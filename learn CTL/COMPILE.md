### Compile standalone executables
Note: makefiles to produce standalone executable are linux specific and will now work under windows

Optimized versions of the software are also available for more high throughput data. 
The follwoing two paragraphs explain how to build a standalone executable version 
using either the C front-end code or the D 2.0 front-end code.

#### (C version)

Clone the repository from github, move into the CTLmappinf folder and run 'make' from a terminal / command line:

    make versionC                                          # Compile the executable
    make static                                            # Compile the static library
    make shared                                            # Compile the shared library

C code can also be compiled also into static or dynamic libraries, when using the 
appropriate makefile commands and can then be used in your own software as external 
dependancy with minimal coupling to your own software.

#### (D 2.0 version)

Prepare your environment by download and installing the DMD 2.0 compiler from 
[www.d-programming-language.org](http://www.d-programming-language.org 
"www.d-programming-language.org"). Run 'make' from a terminal / command line:

    make versionD                                          # Compile the D 2.0 executable

In the future, you can use a provided binary by downloading the approriate one for your 
operating system.
