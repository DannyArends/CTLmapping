#### Start of system configuration section. ####
srcdir = ../Rctl/src
CC = gcc -O -pedantic -w -Wall -Wextra
LINK = ar
DEFS = -DSTANDALONE

CFLAGS = -c $(CDEBUG) -I. -I$(srcdir) $(DEFS)
LDFLAGS = -lm
#### End of system configuration section. ####

SOURCES    = $(srcdir)/vector.c $(srcdir)/matrix.c $(srcdir)/sort.c $(srcdir)/structs.c \
             $(srcdir)/ctlio.c $(srcdir)/correlation.c $(srcdir)/permutation.c $(srcdir)/mapctl.c
OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = mapctl.exe
STATIC     = libctl.a
SHARED     = libctl.so

all: clean $(SOURCES) $(EXECUTABLE)

static: clean $(SOURCES) $(STATIC)

debug: CFLAGS += -g
debug: all

shared: DEFS += -fPIC
shared: clean $(SOURCES) $(SHARED)

$(EXECUTABLE): $(OBJECTS) main.o
	$(CC) $(LDFLAGS) main.o $(OBJECTS) -o $@

$(STATIC): $(OBJECTS)
	$(LINK) rcs $@ $(OBJECTS)

$(SHARED): $(OBJECTS)
	$(CC) -shared -Wl $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o *.so *.a *exe
	rm -rf ../Rctl/src/*o

