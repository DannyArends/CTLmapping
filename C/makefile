#### Start of system configuration section. ####
srcdir = ../Rctl/src
CC = gcc
LINK = ar
DEFS = -DSTANDALONE

CFLAGS = -I/usr/share/R/include -std=gnu99 -g -O3 -pedantic -w -Wall -Wextra -c $(CDEBUG) -I. -I$(srcdir) $(DEFS)
LDFLAGS = -lm -L/usr/lib/R/lib -lR -O3
#### End of system configuration section. ####

SOURCES    = $(srcdir)/vector.c $(srcdir)/matrix.c $(srcdir)/sort.c $(srcdir)/structs.c \
             $(srcdir)/ctlio.c $(srcdir)/correlation.c $(srcdir)/permutation.c \
             $(srcdir)/mapctl.c $(srcdir)/rmapctl.c
OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = mapctl
STATIC     = libctl.a
SHARED     = libctl.so

# target: help - Display possible targets
help:
	egrep "^# target:" [Mm]akefile

all: clean $(SOURCES) $(EXECUTABLE)

# target: static - CTL as static library
static: clean $(SOURCES) $(STATIC)

# target: test - Debug build
debug: CFLAGS += -g
debug: all

# target: test - Test by executing CTL mapping
test: all
	./mapctl

# target: shared - CTL as shared library
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

# target: clean - Cleanup
clean:
	rm -rf *o *.so *.a *exe summary.txt
	rm -rf ../Rctl/src/*o

