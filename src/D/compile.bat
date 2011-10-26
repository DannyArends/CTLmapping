call setpath d
del D_qcl.dll
dmd -ofD_qcl.dll -L/IMPLIB  -O -inline -release D\qcl.d D\dll\dllmain.d D\dll\raux.d D\libs\libload.d D\libs\r.d D\dll\mydll.def phobos.lib -IC:\D\dmd2\windows\bin
del ATA.lib
del D_qcl.obj
