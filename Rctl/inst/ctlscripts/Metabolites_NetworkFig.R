postscript("test.eps", width = 19.0, height = 19.0, paper="special", horizontal=FALSE)
op <- par(mfrow=c(3,3))
for(x in 1:8){  ctl.lineplot(ctls, map_info, x, 0.0001,cex=2, col="darkorange") }
plot(1:10,1:10,t='n',xaxt='n',xlab="",,yaxt='n',ylab="")
for(x in 1:9){
  draw.element(3,10-(x - .5),as.character(x),cex=2, bg.col="lightgray")
  text(6,10-(x - .5),ctl.names(ctls)[x],cex=2)
}
dev.off()
