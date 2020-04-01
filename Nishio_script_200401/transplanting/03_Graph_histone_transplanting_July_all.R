
library(fields)

pdf("figs/FLC_K27_K4_STM_ACT2_180703transplant_cold_all.pdf",width=2.7,height=0.95)
set.panel(1,2,relax=T)
par(oma=c(0,0,0,0))
par(mar=c(1.6,1.6,0.8,0.1))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)
par(xpd=T)

######## K27me3(/STM) cold ##############
y<-read.csv("data/H3K27me3_180703transplant.csv",header=T,sep=",")
y2<-as.matrix(y)
plot(y2[,2],y2[,3],type="o",pch=16,cex=0.5,lwd=1,col=rgb(0,0.7,1),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,2.1))

axis(side=1,at=c(0,10,20,30,40,50),label=F,tcl=-0.2,pos=0)
arrows(-2,0,50,0,length=0)
mtext(side=1,at=c(0,10,20,30,40,50),text=c("0","10","20","30","40","50"),
      line=-0.3)
mtext(side=1,text="Days after transplantation",line=0.1)

axis(side=2,at=seq(0,2,1),tcl=0.2,las=2,pos=-2)
arrows(-2,0,-2,2.1,length=0)
mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.55)
mtext("H3K27me3",side=2,line=0.15)
mtext("To cold (3 Jul.)",side=3,line=-0.1)

par(new=T)
plot(y2[,2],y2[,4],type="o",pch=18,cex=0.5,lwd=1,col=rgb(0,0.2,0.5),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,2.1), lty=1)
par(new=T)
plot(y2[,2],y2[,5],type="o",pch=15,cex=0.5,lwd=1,col=rgb(1,0,0.2),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,2.1), lty=1)
par(new=T)
plot(y2[,2],y2[,6],type="o",pch=17,cex=0.5,lwd=1,col="orange",
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,2.1), lty=1)
par(new=T)
plot(y2[,2],y2[,7],type="o",pch=4,cex=0.5,lwd=1,col=rgb(0.6,0.3,0.6),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,2.1), lty=1)
par(new=T)
plot(y2[,2],y2[,8],type="o",pch=1,cex=0.5,lwd=1,col=rgb(1,0,1),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,2.1), lty=1)

arrows(y2[,2], y2[,3]-y2[,10], y2[,2], y2[,3]+y2[,10], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(0,0.7,1))
arrows(y2[,2], y2[,4]-y2[,11], y2[,2], y2[,4]+y2[,11], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(0,0.2,0.5))
arrows(y2[,2], y2[,5]-y2[,12], y2[,2], y2[,5]+y2[,12], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(1,0,0.2))
arrows(y2[,2], y2[,6]-y2[,13], y2[,2], y2[,6]+y2[,13], length=.01, 
       angle=90, code=3,lwd=0.5, col="orange")
arrows(y2[,2], y2[,7]-y2[,14], y2[,2], y2[,7]+y2[,14], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(0.6,0.3,0.6))
arrows(y2[,2], y2[,8]-y2[,15], y2[,2], y2[,8]+y2[,15], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(1,0,1))

arrows(-2,2.1,50,2.1,length=0)
arrows(50,0,50,2.1,length=0)



######## K4me3(/ACT2) cold ##############
y<-read.csv("data/H3K4me3_180703transplant.csv",header=T,sep=",")
y2<-as.matrix(y)
plot(y2[,2],y2[,3],type="o",pch=16,cex=0.5,lwd=1,col=rgb(0,0.7,1),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,1.05))

axis(side=1,at=c(0,10,20,30,40,50),label=F,tcl=-0.2,pos=0)
arrows(-2,0,50,0,length=0)
mtext(side=1,at=c(0,10,20,30,40,50),text=c("0","10","20","30","40","50"),
      line=-0.3)
mtext(side=1,text="Days after transplantation",line=0.1)

axis(side=2,at=seq(0,1,0.5),tcl=0.2,las=2,pos=-2)
arrows(-2,0,-2,1.05,length=0)
mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.85)
mtext("H3K4me3",side=2,line=0.45)
mtext("To cold (3 Jul.)",side=3,line=-0.1)

par(new=T)
plot(y2[,2],y2[,4],type="o",pch=18,cex=0.5,lwd=1,col=rgb(0,0.2,0.5),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,1.05), lty=1)
par(new=T)
plot(y2[,2],y2[,5],type="o",pch=15,cex=0.5,lwd=1,col=rgb(1,0,0.2),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,1.05), lty=1)
par(new=T)
plot(y2[,2],y2[,6],type="o",pch=17,cex=0.5,lwd=1,col="orange",
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,1.05), lty=1)
par(new=T)
plot(y2[,2],y2[,7],type="o",pch=4,cex=0.5,lwd=1,col=rgb(0.6,0.3,0.6),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,1.05), lty=1)
par(new=T)
plot(y2[,2],y2[,8],type="o",pch=1,cex=0.5,lwd=1,col=rgb(1,0,1),
     axes=F,ann=F,xlim=c(-2,50),ylim=c(0,1.05), lty=1)

arrows(y2[,2], y2[,3]-y2[,10], y2[,2], y2[,3]+y2[,10], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(0,0.7,1))
arrows(y2[,2], y2[,4]-y2[,11], y2[,2], y2[,4]+y2[,11], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(0,0.2,0.5))
arrows(y2[,2], y2[,5]-y2[,12], y2[,2], y2[,5]+y2[,12], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(1,0,0.2))
arrows(y2[,2], y2[,6]-y2[,13], y2[,2], y2[,6]+y2[,13], length=.01, 
       angle=90, code=3,lwd=0.5, col="orange")
arrows(y2[,2], y2[,7]-y2[,14], y2[,2], y2[,7]+y2[,14], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(0.6,0.3,0.6))
arrows(y2[,2], y2[,8]-y2[,15], y2[,2], y2[,8]+y2[,15], length=.01, 
       angle=90, code=3,lwd=0.5, col=rgb(1,0,1))

arrows(-2,1.05,50,1.05,length=0)
arrows(50,0,50,1.05,length=0)


set.panel()
dev.off()
