#### Start of system configuration section. ####
CC = gcc
LINK = ar
#### End of system configuration section. ####

OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = naive
LDFLAGS = -lm

all: clean $(SOURCES) $(EXECUTABLE)

# target: help - Display possible targets
help:
	egrep "^# target:" [Mm]akefile


$(EXECUTABLE): main.o
	$(CC) -O3 main.o $(OBJECTS) $(LDFLAGS) -o $@

# target: clean - Cleanup
clean:
	rm -rf naive

