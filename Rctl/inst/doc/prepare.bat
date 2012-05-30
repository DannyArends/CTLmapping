R CMD Sweave manual.Rnw
latex  manual.tex
bibtex manual
latex  manual.tex
latex  manual.tex
dvipdfm manual.dvi
rm manual.aux
rm manual.tex
rm manual.dvi
rm manual.log
rm manual.bbl
rm manual.blg