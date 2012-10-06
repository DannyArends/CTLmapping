#### Start of system configuration section. ####
srcdir = .
CC = gcc -O
DEFS =
LIBS =

CDEBUG = -g
CFLAGS = $(CDEBUG) -I$(srcdir) $(DEFS)
LDFLAGS = -g
#### End of system configuration section. ####

SOURCES    = main.c vector.c sort.c structs.c ctlio.c correlation.c mapctl.c
OBJECTS    = $(SOURCES:.cpp=.o)
EXECUTABLE = mapctl.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o *exe