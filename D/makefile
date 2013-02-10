#### Start of system configuration section. ####
CC   = rdmd --force --build-only
MAKE = make
MOVE = mv
CSRC = ../C

all: clean buildD

# target buildD: Build the D executable
buildD: libctl.a
	$(CC) -Isrc/ -Llibctl.a src/ctl/mapctl.d
	$(MOVE) src/ctl/mapctl mapctl

# target libctl.a: Create the CTL static library
libctl.a:
	cd $(CSRC); \
	$(MAKE) static; \
	$(MOVE) libctl.a ../D/libctl.a


clean:
	rm -rf *o *.so *.a *.exe *.bin
	rm -rf ../Rctl/src/*o

