#### Start of system configuration section. ####
RCMD = R CMD
MAKE = make
MOVE = mv
CSRC = ./C
DSRC = ./D

all: static shared installR versionC

# target versionC: standalone version from C source code
versionC:
	cd $(CSRC); \
	$(MAKE); \
	$(MOVE) mapctl ../mapctl

# target versionD: standalone version from D source code (QTAB and qtlHD integration)
versionD:
	cd $(DSRC); \
	$(MAKE); \
	$(MOVE) mapctl ../mapctl

# target static: Create the CTL static library
static:
	cd $(CSRC); \
	$(MAKE) static; \
	$(MOVE) libctl.a ../libctl.a

# target shared: Create the CTL shared library
shared:
	cd $(CSRC); \
	$(MAKE) shared; \
	$(MOVE) libctl.so ../libctl.so

# target checkR: Check the R version
checkR: clean
	$(RCMD) check Rctl
	$(MAKE) clean

# target installR: Install the R version
installR: clean
	$(RCMD) INSTALL Rctl
	$(MAKE) clean

# target clean: cleanup
clean:
	rm -f mapctl *a *.so
	rm -rf Rctl/src/*o
	rm -rf Rctl.Rcheck

