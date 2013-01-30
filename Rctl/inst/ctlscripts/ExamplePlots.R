aa <- rnorm(500)
bb <- rnorm(500)
ma <- (runif(500)-0.5)/2
mb <- (runif(500)-0.5)/2

pchaa <- 20
pchbb <- 20

ccex <- 1.0

setwd("~/Github/Articles/CTLpaper/nar/img")
#png("CTL_Example.png",width=1024,height=768)
postscript("CTL_Example.eps",width=20, height=15,paper="special", horizontal=FALSE)
  op <- par(mfrow = c(2,2))
  op <- par(cex=c(2))
  op <- par(mgp=c(0.1, 1, 0))
  op <- par(mai=c(0.2, 0.2, 0.1, 0.1))

  #A) QTL, CTL
  plot(c(-3,3),c(-3,3),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(sort(aa)+ma,sort(bb)+mb, pch=pchaa,col='red')
  points(bb+ma, (0.2*aa)-2.5, pch=pchbb,col='blue')
  legend("topleft", c("A)"),bty='n')

  #B) QTL, CTL
  plot(c(-3,3),c(-3,3),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(aa, bb, pch=pchaa,col='red')
  points(sort(aa)+ma,sort(bb)+mb,pch=pchbb,col="blue")
  legend("topleft", c("B)"),bty='n')

  plot(c(-3,3),c(-3,3),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(bb+ma, (0.2*aa)-1, pch=pchaa,col='red')
  points(bb+ma, (0.2*aa), pch=pchbb,col='blue')
  legend("topleft", c("C)"),bty='n')
  #legend("topright", c("No QTL T1","QTL T2","No CTL T1:T2"),bty='n')


  #legend("topleft", c("No QTL T1","QTL T2","CTL T1:T2"),bty='n')

  plot(c(-3,3),c(-3,3),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(bb, aa, pch=pchaa,col='red')
  points(bb+ma, aa, pch=pchbb,col='blue')
  legend("topleft", c("D)"),bty='n')
  #legend("topright", c("No QTL T1","No QTL T2","No CTL T1:T2"),bty='n')
dev.off()
