#! /bin/sh

echo "* Starting SWEAVE"
rm -vf *.tex
which R
R --version

R CMD BATCH SweaveIt.R
if [ ! -e article.tex ]; then
  cat SweaveIt.Rout 
  exit 2
fi
# latex article.tex
# run twice for references
# latex article.tex 
dvipdfm article.dvi
