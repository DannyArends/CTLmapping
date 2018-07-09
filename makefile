#### Start of system configuration section. ####
RCMD = R CMD
MAKE = make
MOVE = mv
CSRC = ./C
DSRC = ./D
RDOC = ./Rctl/vignettes
RPDF = ./Rctl/inst/doc
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
	cd $(RDOC); \
  $(RCMD) Sweave manual.Rnw; \
  pdflatex manual.tex;\
  rm -rf manual.aux manual.log manual.tex manual.aux.Rcheck
	$(MOVE) $(RDOC)/manual.pdf $(RPDF)
	$(RCMD) check --no-tests Rctl
	$(MAKE) clean

# target: installR - Install the R version
installR: clean
	$(RCMD) INSTALL Rctl
	$(MAKE) clean
 
# target: buildR - Build the R version for CRAN submission
buildR: clean
	$(RCMD) build Rctl
	$(MAKE) clean

# target: clean - Cleanup
clean:
	rm -f mapctl *a *.so q.log
	rm -f C/*o C/summary.txt
	rm -f $(RPDF)/manual.pdf
	rm -rf Rctl/src/*o
	rm -rf Rctl.Rcheck
	cd $(RDOC); \
  rm -f manual.aux manual.log manual.tex
