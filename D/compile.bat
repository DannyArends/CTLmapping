dmd -lib -O -ofarray.lib src/ctl/core/array/matrix.d src/ctl/core/array/search.d src/ctl/core/array/ranges.d
dmd -lib -O -ofstats.lib src/ctl/core/stats/basic.d src/ctl/core/stats/correlation.d src/ctl/core/stats/tolod.d -Isrc/
dmd -lib -O -ofinput.lib src/ctl/io/csv/parse.d src/ctl/io/csv/write.d -Isrc/
dmd -lib -O -ofctl.lib src/ctl/core/ctl/mapping.d src/ctl/core/ctl/permutation.d src/ctl/core/ctl/utils.d -Isrc/
dmd -O src/ctl/mapctl.d array.lib stats.lib input.lib ctl.lib -Isrc/