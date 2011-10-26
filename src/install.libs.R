if(WINDOWS){
  shell(paste(R_PACKAGE_SOURCE,"/src/D/compile.bat",sep=""))
}else{
  shell(paste(R_PACKAGE_SOURCE,"/src/D/compile.sh",sep=""))
}
files <- Sys.glob(paste("*", SHLIB_EXT, sep=''))
libarch <- if (nzchar(R_ARCH)) paste('libs', R_ARCH, sep='') else 'libs'
dest <- file.path(R_PACKAGE_DIR, libarch)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)