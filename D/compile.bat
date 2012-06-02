dmd -lib -O -w -inline -ofarray.lib src/ctl/core/array/matrix.d src/ctl/core/array/search.d src/ctl/core/array/ranges.d
dmd -lib -O -w -inline -ofstats.lib src/ctl/core/stats/basic.d src/ctl/core/stats/correlation.d src/ctl/core/stats/tolod.d -Isrc/
dmd -lib -O -w -inline -ofinput.lib src/ctl/core/analysis.d src/ctl/io/reader.d src/ctl/io/terminal.d src/ctl/io/csv/parse.d src/ctl/io/csv/write.d src/ctl/io/cmdline/parse.d -Isrc/
dmd -lib -O -w -inline -ofctl.lib src/ctl/core/ctl/mapping.d src/ctl/core/ctl/permutation.d src/ctl/core/ctl/utils.d -Isrc/
dmd -lib -O -w -inline -ofqtl.lib src/ctl/core/qtl/LUdecomp.d src/ctl/core/qtl/nrc.d src/ctl/core/qtl/qtl.d src/ctl/core/qtl/regression.d src/ctl/core/qtl/utils.d -Isrc/
dmd -O src/ctl/mapctl.d array.lib stats.lib input.lib ctl.lib qtl.lib -Isrc/
