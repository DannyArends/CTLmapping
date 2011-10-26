call setpath d
del qcl.dll
dmd -ofqcl.dll -L/IMPLIB  -O -inline -release qcl.d dll\dllmain.d dll\raux.d libs\libload.d libs\r.d dll\mydll.def phobos.lib -IC:\D\dmd2\windows\bin
move qcl.dll ../tests
del ATA.lib
del qcl.obj
