doCommand <- function(command){
  if(system(command)==0){
    cat("** DMD compiler found\n")
  }else{
    cat("** WARNING: No DMD compiler found\n")
  }
}

Dfiles <- paste(Sys.glob(c("D/*.d","D/dll/*","D/libs/*")),collapse=" ")

if(WINDOWS){
  doCommand(paste("dmd -ofDcode.dll -L/IMPLIB  -O -inline -release ",Dfiles," phobos.lib -IC:/D/dmd2/windows/bin",sep=""))
}else{
  doCommand(paste("dmd -ofDcode.so ",Dfiles," phobos.lib",sep=""))
}
files <- Sys.glob(paste("*", SHLIB_EXT, sep=''))
libarch <- if (nzchar(R_ARCH)) paste('libs', R_ARCH, sep='') else 'libs'
dest <- file.path(R_PACKAGE_DIR, libarch)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
file.remove(Sys.glob(paste("*", SHLIB_EXT, sep='')))
file.remove(Sys.glob("*.o"))
file.remove(Sys.glob("*.lib"))
file.remove(Sys.glob("*.obj"))