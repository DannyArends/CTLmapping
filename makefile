#### Start of system configuration section. ####
RCMD = R CMD
MAKE = make
MOVE = mv
CSRC = ./C
DSRC = ./D
#### End of system configuration section. ####

all: static shared installR versionC

# target: help - Display possible targets
help:
	egrep "^# target:" [Mm]akefile

# target: versionC - Standalone version from C source code
versionC:
	cd $(CSRC); \
	$(MAKE); \
	$(MOVE) mapctl ../mapctl

# target: versionD - Standalone version from D source code
versionD:
	cd $(DSRC); \
	$(MAKE); \
	$(MOVE) mapctl ../mapctl

# target: static - Create the CTL static library
static:
	cd $(CSRC); \
	$(MAKE) static; \
	$(MOVE) libctl.a ../libctl.a

# target: shared - Create the CTL shared library
shared:
	cd $(CSRC); \
	$(MAKE) shared; \
	$(MOVE) libctl.so ../libctl.so

# target: checkR - Check the R version
checkR: clean
	$(RCMD) check Rctl
	$(MAKE) clean

# target: installR - Install the R version
installR: clean
	$(RCMD) INSTALL Rctl
	$(MAKE) clean

# target: clean - Cleanup
clean:
	rm -f mapctl *a *.so q.log
	rm -rf Rctl/src/*o
	rm -rf Rctl.Rcheck
