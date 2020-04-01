
library(fields)

pdf("figs/FLC_K27_K4_STM_ACT2_A-D.pdf",width=3,height=1.2)
set.panel(1,2,relax=T)
par(oma=c(0,0,0,0))
par(mar=c(0.5,1.3,0.2,0.1))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)

######## K27me3(/STM) ##############
x<-read.csv("data/H3K27me3_H3K4me3_a-d.csv",header=T,sep=",")
y<-read.csv("data/H3K27me3_FLC.csv",header=T,sep=",")
y2<-as.matrix(y)
#z<-read.csv("data/120901-140930temp.csv",header=T,sep=",")
#z2<-as.matrix(z)
#plot(z2[427:577,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,151),ylim=c(-5,35))
#axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=151)
#mtext("Temperature (°C)",side=4,line=0.05)

#par(new=T)
plot(x$days,as.numeric(y2[seq(29,37,2),22]),type="o",pch=18,cex=0.3,lwd=1,
     col=rgb(0,0.2,1),axes=F,ann=F,xlim=c(0,151),ylim=c(0,1.05), lty=1)
par(new=T)
plot(x$days,as.numeric(y2[seq(29,37,2),23]),type="o",pch=15,cex=0.3,lwd=1,
     col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,151),ylim=c(0,1.05), lty=1)

arrows(x$days, as.numeric(y2[seq(29,37,2),22])-as.numeric(y2[seq(29,37,2),31]), 
       x$days, as.numeric(y2[seq(29,37,2),22])+as.numeric(y2[seq(29,37,2),31]),
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.2,1))
arrows(x$days, as.numeric(y2[seq(29,37,2),23])-as.numeric(y2[seq(29,37,2),32]), 
       x$days, as.numeric(y2[seq(29,37,2),23])+as.numeric(y2[seq(29,37,2),32]), 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,0.2))

day<-c(0,30,31,31,28,31)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),sum(day[1:6]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]))
labelx<-c("11","12","1","2","3")
axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
par(font=1)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1,0.5),tcl=0.2,las=2,pos=0)
arrows(0,0,0,1.05,length=0)
#mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.9)
mtext("H3K27me3",side=2,line=0.45)

par(new=T)
plot(x$days,x$K27_a_STM_mean,type="o",pch=16,cex=0.2 ,lty=1,lwd=1,
     axes=F,ann=F,xlim=c(0,151),ylim=c(0,1.05))
par(new=T)
plot(x$days,x$K27_b_STM_mean,type="o",pch=16,cex=0.2,lty=3,lwd=1,
     axes=F,ann=F,xlim=c(0,151),ylim=c(0,1.05))
par(new=T)
plot(x$days,x$K27_c_STM_mean,type="o",pch=16,cex=0.2,lty=5,lwd=1,
     axes=F,ann=F,xlim=c(0,151),ylim=c(0,1.05))
par(new=T)
plot(x$days,x$K27_d_STM_mean,type="o",pch=16,cex=0.2,lty=6,lwd=1,
     axes=F,ann=F,xlim=c(0,151),ylim=c(0,1.05))

arrows(x$days, x$K27_a_STM_mean-x$K27_a_STM_SD, 
       x$days, x$K27_a_STM_mean+x$K27_a_STM_SD, 
       length=.01, angle=90, code=3,lwd=0.5)
arrows(x$days, x$K27_b_STM_mean-x$K27_b_STM_SD, 
       x$days, x$K27_b_STM_mean+x$K27_b_STM_SD, 
       length=.01, angle=90, code=3,lwd=0.5)
arrows(x$days, x$K27_c_STM_mean-x$K27_c_STM_SD, 
       x$days, x$K27_c_STM_mean+x$K27_c_STM_SD, 
       length=.01, angle=90, code=3,lwd=0.5)
arrows(x$days, x$K27_d_STM_mean-x$K27_d_STM_SD, 
       x$days, x$K27_d_STM_mean+x$K27_d_STM_SD, 
       length=.01, angle=90, code=3,lwd=0.5)

arrows(0,1.05,151,1.05,length=0)
arrows(151,0,151,1.05,length=0)


######## K4me3(/ACT2) ##############
x<-read.csv("data/H3K27me3_H3K4me3_A-D.csv",header=T,sep=",")
y<-read.csv("data/H3K4me3_FLC.csv",header=T,sep=",")
y2<-as.matrix(y)
#z<-read.csv("data/120901-140930temp.csv",header=T,sep=",")
#z2<-as.matrix(z)
#plot(z2[427:577,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,151),ylim=c(-5,35))
#axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=151)
#mtext("Temperature (°C)",side=4,line=0.05)

#par(new=T)
plot(x$days,as.numeric(y2[seq(29,37,2),4]),type="o",pch=18,cex=0.3,lwd=1,
     col=rgb(0,0.2,1),axes=F,ann=F,xlim=c(0,151),ylim=c(0,0.9), lty=1)
par(new=T)
plot(x$days,as.numeric(y2[seq(29,37,2),5]),type="o",pch=15,cex=0.3,lwd=1,
     col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,151),ylim=c(0,0.9), lty=1)

arrows(x$days, as.numeric(y2[seq(29,37,2),4])-as.numeric(y2[seq(29,37,2),13]), 
       x$days, as.numeric(y2[seq(29,37,2),4])+as.numeric(y2[seq(29,37,2),13]),
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.2,1))
arrows(x$days, as.numeric(y2[seq(29,37,2),5])-as.numeric(y2[seq(29,37,2),14]), 
       x$days, as.numeric(y2[seq(29,37,2),5])+as.numeric(y2[seq(29,37,2),14]), 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,0.2))

day<-c(0,30,31,31,28,31)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),sum(day[1:6]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]))
labelx<-c("11","12","1","2","3")
axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
par(font=1)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,0.9,0.3),tcl=0.2,las=2,pos=0)
#mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.9)
mtext("H3K4me3",side=2,line=0.45)

par(new=T)
plot(x$days,x$K4_a_ACT2_mean,type="o",pch=16,cex=0.2 ,lty=1,lwd=1,
     axes=F,ann=F,xlim=c(0,151),ylim=c(0,0.9))
par(new=T)
plot(x$days,x$K4_b_ACT2_mean,type="o",pch=16,cex=0.2,lty=3,lwd=1,
     axes=F,ann=F,xlim=c(0,151),ylim=c(0,0.9))
par(new=T)
plot(x$days,x$K4_c_ACT2_mean,type="o",pch=16,cex=0.2,lty=5,lwd=1,
     axes=F,ann=F,xlim=c(0,151),ylim=c(0,0.9))
par(new=T)
plot(x$days,x$K4_d_ACT2_mean,type="o",pch=16,cex=0.2,lty=6,lwd=1,
     axes=F,ann=F,xlim=c(0,151),ylim=c(0,0.9))

arrows(x$days, x$K4_a_ACT2_mean-x$K4_a_ACT2_SD, 
       x$days, x$K4_a_ACT2_mean+x$K4_a_ACT2_SD, 
       length=.01, angle=90, code=3,lwd=0.5)
arrows(x$days, x$K4_b_ACT2_mean-x$K4_b_ACT2_SD, 
       x$days, x$K4_b_ACT2_mean+x$K4_b_ACT2_SD, 
       length=.01, angle=90, code=3,lwd=0.5)
arrows(x$days, x$K4_c_ACT2_mean-x$K4_c_ACT2_SD, 
       x$days, x$K4_c_ACT2_mean+x$K4_c_ACT2_SD, 
       length=.01, angle=90, code=3,lwd=0.5)
arrows(x$days, x$K4_d_ACT2_mean-x$K4_d_ACT2_SD, 
       x$days, x$K4_d_ACT2_mean+x$K4_d_ACT2_SD, 
       length=.01, angle=90, code=3,lwd=0.5)

arrows(0,0.9,151,0.9,length=0)
arrows(151,0,151,0.9,length=0)

set.panel()
dev.off()





pdf("figs/FLC_K27_K4_STM_ACT2_A-D_legend.pdf",width=0.35,height=0.8)
par(mar=c(0.1,0,0.1,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)

plot(NULL,xlim=c(0,200),ylim=c(0,1),axes=F,ann=F)

legend(-50,1.2,
       legend="II",lty=1,lwd=1,col=rgb(0,0.2,1),bty="n",x.intersp=0.1,seg.len=1)

legend(-50,1.05,
       legend="A",lty=1,lwd=1,bty="n",x.intersp=0.1,seg.len=1)
legend(-50,0.9,
       legend="B",lty=3,lwd=1,bty="n",x.intersp=0.1,seg.len=1)
legend(-50,0.75,
       legend="C",lty=5,lwd=1,bty="n",x.intersp=0.1,seg.len=1)
legend(-50,0.6,
       legend="D",lty=6,lwd=1,bty="n",x.intersp=0.1,seg.len=1)

legend(-50,0.45,
       legend="III",lty=1,lwd=1,col=rgb(1,0,0.2),bty="n",x.intersp=0.1,seg.len=1)

dev.off()
