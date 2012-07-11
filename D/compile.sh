dmd -lib -O -w -inline -ofarray.a src/ctl/core/array/matrix.d src/ctl/core/array/search.d src/ctl/core/array/ranges.d -I~/github/phobos -I~/github/druntime/import
echo "array"
dmd -lib -O -w -inline -ofstats.a src/ctl/core/stats/basic.d src/ctl/core/stats/correlation.d src/ctl/core/stats/tolod.d -Isrc/ -I~/github/phobos -I~/github/druntime/import
echo "stats"
dmd -lib -O -w -inline -ofinput.a src/ctl/core/analysis.d src/ctl/io/reader.d src/ctl/io/terminal.d src/ctl/io/csv/parse.d src/ctl/io/csv/write.d src/ctl/io/cmdline/parse.d -Isrc/ -I~/github/phobos -I~/github/druntime/import
echo "input"
dmd -lib -O -w -inline -ofctl.a src/ctl/core/ctl/mapping.d src/ctl/core/ctl/permutation.d src/ctl/core/ctl/utils.d -Isrc/ -I~/github/phobos -I~/github/druntime/import
echo "mapping"
dmd -lib -O -w -inline -ofqtl.a src/ctl/core/qtl/LUdecomp.d src/ctl/core/qtl/nrc.d src/ctl/core/qtl/qtl.d src/ctl/core/qtl/regression.d src/ctl/core/qtl/utils.d -Isrc/ -I~/github/phobos -I~/github/druntime/import
echo "qtl"
dmd -O src/ctl/mapctl.d array.a stats.a input.a ctl.a qtl.a -Isrc/ -I~/github/phobos -I~/github/druntime/import -L~/github/phobos/generated/linux/release/32
echo "done"
