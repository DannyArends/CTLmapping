gcc -c vector.c
gcc -c structs.c
gcc -c ctlio.c
gcc -c correlation.c
gcc -c sort.c
gcc -c mapctl.c
gcc -W -Wall main.c vector.o structs.o sort.o ctlio.o correlation.o mapctl.o -o mapctl.exe
