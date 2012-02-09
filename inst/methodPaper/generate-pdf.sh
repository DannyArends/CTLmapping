#! /bin/sh
#
# Create PDF - run this script several times to get the citations
# correct.

bibtex document
# perl nat2jour.pl -maxauth 2 BioGem
latex document.tex
if [ "$1" = "--pdf" ] ; then
  dvipdf document.dvi document.pdf
else
  echo "Skipped generating PDF!"
fi
