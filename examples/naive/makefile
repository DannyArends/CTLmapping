#### Start of system configuration section. ####
CC = gcc
LINK = ar
#### End of system configuration section. ####

OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = naive

all: clean $(SOURCES) $(EXECUTABLE)

# target: help - Display possible targets
help:
	egrep "^# target:" [Mm]akefile


$(EXECUTABLE): main.o
	$(CC) $(LDFLAGS) main.o $(OBJECTS) -o $@

# target: clean - Cleanup
clean:
	rm -rf naive

