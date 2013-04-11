aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2

pchaa <- 20
pchbb <- 20

ccex <- 1.0

setwd("~/Github/Articles/CTLpaper/nar/img")
#png("CTL_Example.png",width=1024,height=768)
postscript("CTL_Example_11.eps",width=20, height=20,paper="special", horizontal=FALSE)
  op <- par(mfrow = c(3,3))
  op <- par(cex=c(2))
  op <- par(mgp=c(0.1, 0, 1))
  op <- par(mai=c(0.1, 0.8, 0.8, 0.1))

aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2

  #A) QTL QTL, CTL
  plot(c(-3.5,3.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(sort(aa)+0.5, bb-1, pch=pchaa,col='black')
  points(sort(aa)+ma-1,sort(bb)+mb+.5,pch=pchbb,col="orange")
  legend("topleft", c("A)"),bty='n')
  mtext("CTL", side=3, line=1, las = 1, cex = 2, font=2)
  mtext("QTL in both traits", side=2, line=1, cex = 2, font=2)
#  op <- par(mai=c(0.1, 0.1, 1, 0.1))

aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2

  plot(c(-3.5,3.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(bb+ma+0.7, (0.2*aa)-0.7+(aa/2), pch=pchaa,col='black')
aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2
  points(bb+ma-0.7, (0.2*aa)+(aa/2)+0.7, pch=pchbb,col='orange')
  legend("topleft", c("D)"),bty='n')
  mtext("No CTL", side=3, line=1, las = 1, cex = 2, font=2)
#aa <- rnorm(1000)
#bb <- rnorm(1000)
#ma <- (runif(1000)-0.5)/2
#mb <- (runif(1000)-0.5)/2

  #B) QTL, CTL
#  plot(c(-3.5,3.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
#  points(sort(aa)-1, bb, pch=pchaa,col='black')
#  points(sort(aa)+ma+0.5,sort(bb)+mb,pch=pchbb,col="orange")
#  legend("topleft", c("B)"),bty='n')

  #A) QTL, CTL - Genotype BB OFF
plot(c(-3.5,4.5),c(-3.5,4.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
points(sort(aa)+ma+1.7,sort(bb)+mb+1.7, pch=pchaa,col='orange')
points(abs(0.1*aa)-2.3,abs(0.1*bb)-2.3, pch=pchbb,col='black')
legend("topleft", c("E)"),bty='n')
mtext("False CTL", side=3, line=1, las = 1, cex = 2, font=2)

aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2

#op <- par(mai=c(0.1, 1, 0.1, 0.1))

  #B) QTL, CTL
  plot(c(-3.5,3.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(sort(aa), bb-1, pch=pchaa,col='black')
  points(sort(aa)+ma,sort(bb)+mb+.5,pch=pchbb,col="orange")
  legend("topleft", c("B)"),bty='n')
  mtext("QTL in one trait", side=2, line=1, cex = 2, font=2)
#op <- par(mai=c(0.1, 0.1, 0.1, 0.1))

#aa <- rnorm(1000)
#bb <- rnorm(1000)
#ma <- (runif(1000)-0.5)/2
#mb <- (runif(1000)-0.5)/2

#  plot(c(-3.5,3.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
#  points(bb+ma, (0.2*aa)-0.7+(aa/2), pch=pchaa,col='black')
#aa <- rnorm(1000)
#bb <- rnorm(1000)
#ma <- (runif(1000)-0.5)/2
#mb <- (runif(1000)-0.5)/2
#  points(bb+ma, (0.2*aa)+(aa/2)+0.7, pch=pchbb,col='orange')
#  legend("topleft", c("D)"),bty='n')
  #legend("topright", c("No QTL T1","QTL T2","No CTL T1:T2"),bty='n')

  plot(c(-3.5,3.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(bb+ma+0.7, (0.2*aa)+(aa/2), pch=pchaa,col='black')
aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2
  points(bb+ma-0.7, (0.2*aa)+(aa/2), pch=pchbb,col='orange')
  legend("topleft", c("D)"),bty='n')
  #legend("topright", c("No QTL T1","QTL T2","No CTL T1:T2"),bty='n')


aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2

  #A) QTL, CTL - Genotype AA OFF
  plot(c(-3.5,3.5),c(-3.5,4.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(sort(aa)+ma,sort(bb)+mb+1.3, pch=pchaa,col='orange')
  points(bb+ma, rep(-2, 1000), pch=pchbb,col='black')
  legend("topleft", c("E)"),bty='n')


#op <- par(mai=c(0.1, 1, 0.1, 0.1))

aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2

  #B) QTL, CTL - Genotype AA ON
  plot(c(-3.5,3.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(aa, bb, pch=pchaa,col='black')
  points(sort(aa)+ma,sort(bb)+mb,pch=pchbb,col="orange")
  legend("topleft", c("C)"),bty='n')
  mtext("No QTL", side=2, line=1, cex = 2, font=2)
#op <- par(mai=c(0.1, 0.1, 0.1, 0.1))

  plot(c(-3.5,3.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  points(bb, aa, pch=pchaa,col='black')

aa <- rnorm(1000)
bb <- rnorm(1000)
ma <- (runif(1000)-0.5)/2
mb <- (runif(1000)-0.5)/2

  points(bb+ma, aa, pch=pchbb,col='orange')
  points(bb[1:50], aa[1:50], pch=pchaa,col='black')
  legend("topleft", c("D)"),bty='n')


#plot(c(-3.5,4.5),c(-3.5,4.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')

  #A) QTL, CTL - Genotype BB OFF
#  plot(c(-3.5, 4.5),c(-3.5,3.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
#  points(sort(aa)+ma+1.3,sort(bb)+mb+.3, pch=pchaa,col='black')
#  points(rep(-2, 1000),bb+ma, pch=pchbb,col='orange')
#  legend("topleft", c("E)"),bty='n')

  #A) QTL, CTL - Genotype AA an BB OFF
#  plot(c(-3.5,4.5),c(-3.5,4.5),t='n', xlab="", ylab="", main="",xaxt='n',yaxt='n')
  #points((0.1*aa)-2.5,bb+ma+1, pch=pchaa,col='orange')
  #points(bb+ma+1, rep(-2, 1000), pch=pchbb,col='black')
  #legend("topleft", c("E)"),bty='n')



  #legend("topleft", c("No QTL T1","QTL T2","CTL T1:T2"),bty='n')

#legend("topright", c("No QTL T1","No QTL T2","No CTL T1:T2"),bty='n')
dev.off()
