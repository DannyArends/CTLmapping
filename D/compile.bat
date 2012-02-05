dmd -lib -O -ofarray.lib src/qcl/core/array/matrix.d src/qcl/core/array/search.d src/qcl/core/array/ranges.d
dmd -lib -O -ofstats.lib src/qcl/core/stats/basic.d src/qcl/core/stats/correlation.d src/qcl/core/stats/tolod.d -Isrc/
dmd -lib -O -ofinput.lib src/qcl/io/csv/parse.d src/qcl/io/csv/write.d -Isrc/
dmd -lib -O -ofqcl.lib src/qcl/core/qcl/singleqcl.d src/qcl/core/qcl/permuteqcl.d src/qcl/core/qcl/utils.d -Isrc/
dmd -O src/qcl/mapqcl.d array.lib stats.lib input.lib qcl.lib -Isrc/