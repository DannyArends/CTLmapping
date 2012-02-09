set QTLHD=E:/github/qtlHD/src/D
dmd -lib -O -ofqtab.lib src/ctl/io/qtab/wrapper.d %QTLHD%/qtl/core/primitives.d %QTLHD%/qtl/core/chromosome.d %QTLHD%/qtl/core/genotype.d %QTLHD%/qtl/core/phenotype.d %QTLHD%/qtl/plugins/qtab/read_qtab.d -Isrc/
dmd -lib -O -ofinput.lib src/ctl/io/reader.d src/ctl/io/csv/parse.d src/ctl/io/csv/write.d src/ctl/io/cmdline/parse.d qtab.lib -Isrc/ -I%QTLHD%/ -version=QTAB
dmd -O src/ctl/mapctl.d array.lib stats.lib input.lib ctl.lib -Isrc/