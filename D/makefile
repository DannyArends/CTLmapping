#### Start of system configuration section. ####
CC   = rdmd --force --build-only
MAKE = make
MOVE = mv
CSRC = ../C
#### End of system configuration section. ####

all: clean buildD

# target: help - Display possible targets
help:
	egrep "^# target:" [Mm]akefile

# target: buildD - Build the D executable (DEFAULT)
# We need to load the R shared library so we link with -L-ldl -L-lR
buildD: libctl
	$(CC) -Isrc/ -Llibctl.a -L-ldl -L-lR src/ctl/mapctl.d
	$(MOVE) src/ctl/mapctl mapctl

# target: libctl - Create the CTL static library
libctl:
	cd $(CSRC); \
	$(MAKE) static; \
	$(MOVE) libctl.a ../D/libctl.a

# target: test - Test by executing CTL mapping
test: all
	./mapctl

# target: clean - Cleanup
clean:
	rm -rf *o *.so *.a *.exe *.bin
	rm -rf ../Rctl/src/*o

