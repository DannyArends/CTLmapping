#### Start of system configuration section. ####
srcdir = ../Rctl/src
CC = gcc
LINK = ar
DEFS = -DSTANDALONE
RINCLUDE = "${R_HOME}/include/"

ifeq ($(OS),Windows_NT)
	RLIB = "${R_HOME}/bin/x64/"
else
	RLIB = "/usr/lib/"
endif

CFLAGS = -I $(RINCLUDE) -std=gnu99 -O3 -pedantic -w -Wall -Wextra -c $(CDEBUG) -I. -I$(srcdir) $(DEFS)
LDFLAGS = -lm -L $(RLIB) -lR -O3
#### End of system configuration section. ####

SOURCES    = $(srcdir)/vector.c $(srcdir)/matrix.c $(srcdir)/sort.c $(srcdir)/structs.c \
             $(srcdir)/ctlio.c $(srcdir)/correlation.c $(srcdir)/permutation.c \
             $(srcdir)/mapctl.c $(srcdir)/rmapctl.c
OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = mapctl
STATIC     = libctl.a
SHARED     = libctl.so

all: clean $(SOURCES) $(EXECUTABLE)

# target: help - Display possible targets
help:
	egrep "^# target:" [Mm]akefile

# target: static - CTL as static library
static: clean $(SOURCES) $(STATIC)

# target: test - Debug build
debug: CFLAGS += -g -DTEST
debug: LDFLAGS += -g -DTEST
debug: all

# target: test - Test by executing CTL mapping
test: debug
#	./mapctl -h
#	./mapctl

# target: shared - CTL as shared library
shared: DEFS += -fPIC
shared: clean $(SOURCES) $(SHARED)

$(EXECUTABLE): $(OBJECTS) main.o
	$(CC)  main.o $(OBJECTS) $(LDFLAGS) -o $@

$(STATIC): $(OBJECTS)
	$(LINK) rcs $@ $(OBJECTS)

$(SHARED): $(OBJECTS)
	$(CC) -shared  $(OBJECTS) $(LDFLAGS)  -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

# target: clean - Cleanup
clean:
	rm -rf *o *.so *.a *exe summary.txt
	rm -rf ../Rctl/src/*o

