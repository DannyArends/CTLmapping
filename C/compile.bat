gcc -c structs.c
gcc -c ctlio.c
gcc -c correlation.c
gcc -c mapctl.c
gcc -W -Wall main.c structs.o ctlio.o correlation.o mapctl.o
