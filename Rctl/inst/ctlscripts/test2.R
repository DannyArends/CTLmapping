
plot.CTLscan3 <- function(x, map_info, todo = c(4,6,2),xaxt='n',xlab=""){
  p <- rep(0,ncol(x$l))
  i <- 1;
  summarized <- apply(x$l,2,sum)
  xxx <- NULL
  mycolors <- topo.colors(nrow(x$l[todo,]))
  plot(c(0, chr_total_length(map_info,gap=0)),c(0, 100), type='n',xlab=xlab,ylab="LOD score",xaxt=xaxt)
  chredge <- 0
  colo <- "white"
  for(f in unique(map_info[,1])){ nchredge <- chredge + chrlength(map_info, f); 
    rect(chredge,-10,nchredge,110, col=colo,border=NA)
    chredge <- nchredge
    if(colo == "gray"){ colo <- "white"; }else{ colo <- "gray"; }
  }
  apply(x$l,1,
    function(d){
     for(idx in 1:length(d)){
     #   rect(m_loc(map_info,idx)-1,p[idx],m_loc(map_info,idx)+1,p[idx]+d[idx],col=mycolors[i],lwd=0,lty=0)
      }
      p <<- p + d
      i <<- i + 1
    }
  )
  loc <- NULL
  for(y in 1:nrow(map_info)){loc <- c(loc,m_loc(map_info,y,gap=0))}

  for(y in unique(map_info[,1])){
    onchr <- which(map_info[,1]==y)
    points(loc[onchr],summarized[onchr],type='l',lwd=2,lty=1)
    points(loc[onchr],as.numeric(x$qtl)[onchr],type='l',lwd=2,lty=2,col="red")
  }
  box()
}

#[c(4,2,6,1,3,5,7,8,9),]

#pathway
postscript("CTL_heatmaps.eps",width=30, height=10,paper="special")
op <- par(mfrow=c(3,2))
op <- par(mgp=c(1, 0.2, 0))
op <- par(mai=c(0.5, 2.0, 0.1, 0.1))
op <- par(cex=c(1.5))
#MT4
plot.CTLscan3(ctls[[6]],map_info)
legend("topleft",c("QTL profile MT4","CTL profile MT4"),col=c("red","black"),lwd=2,lty=c(2,1),bty="n")
image(c(1:69),c(1:9),t(ctls[[6]]$l[c(4,2,6,1,3,5,7,8,9),]),breaks=c(0,3,5,10,25,1000000),col=c("white","lightgray","gray","darkgray","black"),yaxt='n',xlab="",ylab="",xaxt='n')
axis(2,gsub(".Mean","  ",colnames(metabolites))[c(4,2,6,1,3,5,7,8,9)],at=c(1:9),las=2,cex.axis=0.7)
for(x in 1:3){rect(0,0,69.5,x+.5,border='blue',lwd=0.1)}
for(x in 4:6){rect(0,0,69.5,x+.5,border='red',lwd=0.1)}
for(x in 7:8){rect(0,0,69.5,x+.5,border='green',lwd=0.1)}
legend("topleft",c("3 < LOD < 5","5 < LOD < 10","10 < LOD < 25","LOD 25+"),col=c("lightgray","gray","darkgray","black"),lwd=10,lty=c(1,1),bg="white")
box()
#MSO4
plot.CTLscan3(ctls[[2]],map_info)
legend("topleft",c("QTL profile MSO4","CTL profile MSO4"),col=c("red","black"),lwd=2,lty=c(2,1),bty="n")
image(c(1:69),c(1:9),t(ctls[[2]]$l[c(4,2,6,1,3,5,7,8,9),]),breaks=c(0,3,5,10,25,1000000),col=c("white","lightgray","gray","darkgray","black"),yaxt='n',xlab="",ylab="",xaxt='n')
axis(2,gsub(".Mean","  ",colnames(metabolites))[c(4,2,6,1,3,5,7,8,9)],at=c(1:9),las=2,cex.axis=0.7)
for(x in 1:3){rect(0,0,69.5,x+.5,border='blue',lwd=0.1)}
for(x in 4:6){rect(0,0,69.5,x+.5,border='red',lwd=0.1)}
for(x in 7:8){rect(0,0,69.5,x+.5,border='green',lwd=0.1)}
box()
#But-3-enyl
plot.CTLscan3(ctls[[4]],map_info,xaxt='t',xlab="Location (cM)")
legend("topleft",c("QTL profile But-3-enyl","CTL profile But-3-enyl"),col=c("red","black"),lwd=2,lty=c(2,1),bty="n")
image(c(1:69),c(1:9),t(ctls[[4]]$l[c(4,2,6,1,3,5,7,8,9),]),breaks=c(0,3,5,10,25,1000000),col=c("white","lightgray","gray","darkgray","black"),yaxt='n',xlab="Marker",ylab="")
axis(2,gsub(".Mean","  ",colnames(metabolites))[c(4,2,6,1,3,5,7,8,9)],at=c(1:9),las=2,cex.axis=0.7)
for(x in 1:3){rect(0,0,69.5,x+.5,border='blue',lwd=0.1)}
for(x in 4:6){rect(0,0,69.5,x+.5,border='red',lwd=0.1)}
for(x in 7:8){rect(0,0,69.5,x+.5,border='green',lwd=0.1)}
box()
dev.off()



#pathway
op <- par(mfrow=c(3,2))
op <- par(mgp=c(2.5, 1, 0))
op <- par(mai=c(0.5, 0.5, 0.1, 0.1))
#MT3
plot.CTLscan3(ctls[[5]],map_info)
image(c(1:69),c(1:3),t(ctls[[5]]$l[c(1,3,5),]),breaks=c(0,3,10,1000),col=c("white","gray","black"),yaxt='n',xlab="Marker",ylab="")
axis(2,c("OPH3","Allyl","MT3"),at=c(1,2,3),las=2,cex.axis=0.7)
box()
#MSO4
plot.CTLscan3(ctls[[3]],map_info)
image(c(1:69),c(1:9),t(ctls[[3]]$l[c(1,3,5),]),breaks=c(0,3,10,1000),col=c("white","gray","black"),yaxt='n',xlab="Marker",ylab="")
axis(2,ctl.names(ctls),at=c(1:9),las=2,cex.axis=0.7)
box()
#But-3-enyl
plot.CTLscan3(ctls[[1]],map_info)
image(c(1:69),c(1:3),t(ctls[[1]]$l[c(1,3,5),]),breaks=c(0,3,10,1000),col=c("white","gray","black"),yaxt='n',xlab="Marker",ylab="")
axis(2,c("OPH3","Allyl","MT3"),at=c(1,2,3),las=2,cex.axis=0.7)
box()



#pathway
op <- par(mfrow=c(4,2))
op <- par(mgp=c(2.5, 1, 0))
op <- par(mai=c(0.5, 0.5, 0.1, 0.1))
#MT3
plot.CTLscan3(ctls[[5]],map_info)
image(c(1:69),c(1:4),t(ctls[[5]]$l[c(8,9,6,5),]),breaks=c(0,3,10,1000),col=c("white","gray","black"),yaxt='n',xlab="Marker",ylab="")
axis(2,c("MT8","MT7","MT4","MT3"),at=c(1,2,3,4),las=2,cex.axis=0.7)
box()
#MT4
plot.CTLscan3(ctls[[6]],map_info)
image(c(1:69),c(1:4),t(ctls[[6]]$l[c(8,9,6,5),]),breaks=c(0,3,10,1000),col=c("white","gray","black"),yaxt='n',xlab="Marker",ylab="")
axis(2,c("MT8","MT7","MT4","MT3"),at=c(1,2,3,4),las=2,cex.axis=0.7)
box()
#MT7
plot.CTLscan3(ctls[[9]],map_info)
image(c(1:69),c(1:4),t(ctls[[9]]$l[c(8,9,6,5),]),breaks=c(0,3,10,1000),col=c("white","gray","black"),yaxt='n',xlab="Marker",ylab="")
axis(2,c("MT8","MT7","MT4","MT3"),at=c(1,2,3,4),las=2,cex.axis=0.7)
box()
#MT8
plot.CTLscan3(ctls[[8]],map_info)
image(c(1:69),c(1:4),t(ctls[[8]]$l[c(8,9,6,5),]),breaks=c(0,3,10,1000),col=c("white","gray","black"),yaxt='n',xlab="Marker",ylab="")
axis(2,c("MT8","MT7","MT4","MT3"),at=c(1,2,3,4),las=2,cex.axis=0.7)
box()


