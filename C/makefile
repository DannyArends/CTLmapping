#### Start of system configuration section. ####
srcdir = ../Rctl/src
CC = gcc -O -pedantic -w -Wall -Wextra
DEFS =
LIBS =

CDEBUG = -g
CFLAGS = $(CDEBUG) -I. -I$(srcdir) $(DEFS)
LDFLAGS = -g -I. -I$(srcdir) -lm -DSTANDALONE
#### End of system configuration section. ####

SOURCES    = main.c $(srcdir)/vector.c $(srcdir)/matrix.c $(srcdir)/sort.c $(srcdir)/structs.c \
                    $(srcdir)/ctlio.c $(srcdir)/correlation.c $(srcdir)/permutation.c $(srcdir)/mapctl.c
OBJECTS    = $(SOURCES:.cpp=.o)
EXECUTABLE = mapctl.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o *exe
