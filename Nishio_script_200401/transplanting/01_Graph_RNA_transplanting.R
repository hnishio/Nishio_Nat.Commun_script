
library(fields)

dir.create("figs")

pdf("figs/FLC_RNA_transplant.pdf",width=3.3,height=0.9)
par(mar=c(1.2,2.1,0.4,4.5))
par(cex=1)
par(ps=6)
par(mgp=c(0,0.1,0))
par(xpd=T)

RNA<-read.csv("data/RNA_natural_2018.csv",header=T,sep=",")
plot(RNA$date+13,RNA$Ave.log10.FLC.ACT2,type="o",pch=16,cex=0.3,lwd=1,
     col="gray",axes=F,ann=F,xlim=c(0,212),ylim=c(-3,1))
arrows(RNA$date+13, RNA$Ave.log10.FLC.ACT2-RNA$Std.log10.FLC.ACT2, 
       RNA$date+13, RNA$Ave.log10.FLC.ACT2+RNA$Std.log10.FLC.ACT2, 
       length=.01, angle=90, code=3,lwd=0.5, col="gray")

par(new=T)
RNA2<-read.csv("data/RNA_180326transplant.csv",header=T,sep=",")
plot(RNA2$date+54,RNA2$Ave.log10.FLC.ACT2.W,type="o",pch=16,cex=0.3,lwd=1,
     col=rgb(1,0.65,0),axes=F,ann=F,xlim=c(0,212),ylim=c(-3,1))
arrows(RNA2$date+54, RNA2$Ave.log10.FLC.ACT2.W-RNA2$Std.log10.FLC.ACT2.W, 
       RNA2$date+54, RNA2$Ave.log10.FLC.ACT2.W+RNA2$Std.log10.FLC.ACT2.W, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0.65,0))

par(new=T)
plot(RNA2$date+54,RNA2$Ave.log10.FLC.ACT2.C,type="o",pch=16,cex=0.3,lwd=1,
     col="skyblue",axes=F,ann=F,xlim=c(0,212),ylim=c(-3,1))
arrows(RNA2$date+54, RNA2$Ave.log10.FLC.ACT2.C-RNA2$Std.log10.FLC.ACT2.C, 
       RNA2$date+54, RNA2$Ave.log10.FLC.ACT2.C+RNA2$Std.log10.FLC.ACT2.C, 
       length=.01, angle=90, code=3,lwd=0.5, col="skyblue")

par(new=T)
RNA3<-read.csv("data/RNA_180410transplant.csv",header=T,sep=",")
plot(RNA3$date+69,RNA3$Ave.log10.FLC.ACT2.C,type="o",pch=16,cex=0.3,lwd=1,
     col="blue",axes=F,ann=F,xlim=c(0,212),ylim=c(-3,1))
arrows(RNA3$date+69, RNA3$Ave.log10.FLC.ACT2.C-RNA3$Std.log10.FLC.ACT2.C, 
       RNA3$date+69, RNA3$Ave.log10.FLC.ACT2.C+RNA3$Std.log10.FLC.ACT2.C, 
       length=.01, angle=90, code=3,lwd=0.5, col="blue")

par(new=T)
RNA4<-read.csv("data/RNA_180424transplant.csv",header=T,sep=",")
plot(RNA4$date+83,RNA4$Ave.log10.FLC.ACT2.C,type="o",pch=16,cex=0.3,lwd=1,
     col="darkblue",axes=F,ann=F,xlim=c(0,212),ylim=c(-3,1))
arrows(RNA4$date+83, RNA4$Ave.log10.FLC.ACT2.C-RNA4$Std.log10.FLC.ACT2.C, 
       RNA4$date+83, RNA4$Ave.log10.FLC.ACT2.C+RNA4$Std.log10.FLC.ACT2.C, 
       length=.01, angle=90, code=3,lwd=0.5, col="darkblue")

par(new=T)
RNA5<-read.csv("data/RNA_180703transplant.csv",header=T,sep=",")
plot(RNA5$date+153,RNA5$Ave.log10.FLC.ACT2.C,type="o",pch=16,cex=0.3,lwd=1,
     col="palegreen4",axes=F,ann=F,xlim=c(0,212),ylim=c(-3,1))
arrows(RNA5$date+153, RNA5$Ave.log10.FLC.ACT2.C-RNA5$Std.log10.FLC.ACT2.C, 
       RNA5$date+153, RNA5$Ave.log10.FLC.ACT2.C+RNA5$Std.log10.FLC.ACT2.C, 
       length=.01, angle=90, code=3,lwd=0.5, col="palegreen4")

day<-c(0,28,31,30,31,30,31,31)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),sum(day[1:6]),sum(day[1:7]),sum(day[1:8]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]),mean(place[6:7]),mean(place[7:8]))
labelx<-c("2","3","4","5","6","7","8")
axis(side=1,at=place,label=F,tcl=-0.2,pos=-3)
mtext(side=1,at=lit,labelx,line=-0.45)
mtext("Month",side=1,line=-0.05)

axis(side=2,at=seq(-3,1,1),tcl=0.2,las=2,pos=0,labels=F)
arrows(0,-3,0,1,length=0)
mtext(side=2,at=seq(-3,1,1),text=c("0.001","0.01","0.1","1","10"),las=1,line=-0.3)
mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.85)
mtext("mRNA",side=2,line=0.5)

legend(212,2.2,legend="Natural",
       col="gray",lty=1,pch=16,lwd=1,pt.cex=0.3,bty="n",x.intersp=0.2,seg.len=0.7)
legend(212,1.5,legend="To warm (26 Mar.)",
       col=rgb(1,0.65,0),lty=1,pch=16,lwd=1,pt.cex=0.3,bty="n",x.intersp=0.2,seg.len=0.7)
legend(212,0.8,legend="To cold (26 Mar.)",
       col="skyblue",lty=1,pch=16,lwd=1,pt.cex=0.3,bty="n",x.intersp=0.2,seg.len=0.7)
legend(212,0.1,legend="To cold (10 Apr.)",
       col="blue",lty=1,pch=16,lwd=1,pt.cex=0.3,bty="n",x.intersp=0.2,seg.len=0.7)
legend(212,-0.6,legend="To cold (24 Apr.)",
       col="darkblue",lty=1,pch=16,lwd=1,pt.cex=0.3,bty="n",x.intersp=0.2,seg.len=0.7)
legend(212,-1.3,legend="To cold (3 Jul.)",
       col="palegreen4",lty=1,pch=16,lwd=1,pt.cex=0.3,bty="n",x.intersp=0.2,seg.len=0.7)

arrows(0,1,212,1,length=0)
arrows(212,-3,212,1,length=0)
dev.off()
