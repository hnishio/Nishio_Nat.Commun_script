
library(fields)

dir.creat("figs")

day<-c(0,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),sum(day[1:11]),sum(day[1:12]),sum(day[1:13]),sum(day[1:14]),sum(day[1:15]),sum(day[1:16]),sum(day[1:17]),sum(day[1:18]),sum(day[1:19]),sum(day[1:20]),sum(day[1:21]),sum(day[1:22]),sum(day[1:23]),sum(day[1:24]),sum(day[1:25]),sum(day[1:26]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]),mean(place[13:14]),mean(place[14:15]),mean(place[15:16]),mean(place[16:17]),mean(place[17:18]),mean(place[18:19]),mean(place[19:20]),mean(place[20:21]),mean(place[21:22]),mean(place[22:23]),mean(place[23:24]),mean(place[24:25]),mean(place[25:26]))
#labelx<-c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep")
labelx<-c("9","10","11","12","1","2","3","4","5","6","7","8","9","10","11","12","1","2","3","4","5","6","7","8","9")


pdf("figs/FLC_K27_K4_RNA_STM_ACT2.pdf",width=3.7,height=4.33)
set.panel(5,1,relax=T)
par(oma=c(0.8,0,0.6,0))
par(mar=c(0.55,2.1,0,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)


###### RNA FLC/ACT2 #######################
z<-read.csv("data/120901-140930temp.csv",header=T,sep=",")
z2<-as.matrix(z)

plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

par(new=T)
x<-read.csv("data/RNA_FLC.csv",header=T,sep=",")
x2<-as.matrix(x)
plot(x2[,"days"],x2[,"FLC.ACT2.log10.mean"],type="o",pch=16,cex=0.5,
     lwd=1,col=rgb(0,0,0),axes=F,ann=F,xlim=c(0,760),ylim=c(-4,2))
arrows(x2[,"days"], x2[,"FLC.ACT2.log10.mean"]-x2[,"FLC.ACT2.log10.SD"], 
       x2[,"days"], x2[,"FLC.ACT2.log10.mean"]+x2[,"FLC.ACT2.log10.SD"], 
       length=.01, angle=90, code=3,lwd=0.5)

axis(side=1,at=place,label=F,tcl=-0.2,pos=-4)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(-4,2,2),lab=c("0.0001","0.01","1","100"),tcl=0.2,las=2,pos=0)
mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.65)
mtext("mRNA",side=2,line=0.3)

arrows(0,2,760,2,length=0)


######## K4me3(/ACT2) ##############
plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

x<-read.csv("data/H3K4me3_FLC.csv",header=T,sep=",")

par(new=T)
plot(x[,"days"],x[,"K4.I.ACT2.mean"],type="o",pch=16,cex=0.5,lwd=1,
     col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5))
par(new=T)
plot(x[,"days"],x[,"K4.II.ACT2.mean"],type="o",pch=18,cex=0.5,lwd=1,
     col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5))
par(new=T)
plot(x[,"days"],x[,"K4.III.ACT2.mean"],type="o",pch=15,cex=0.5,lwd=1,
     col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5))
par(new=T)
plot(x[,"days"],x[,"K4.IV.ACT2.mean"],type="o",pch=17,cex=0.5,lwd=1,
     col="orange",axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5))
par(new=T)
plot(x[,"days"],x[,"K4.V.ACT2.mean"],type="o",pch=4,cex=0.5,lwd=1,
     col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5))

axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1.5,0.5),tcl=0.2,las=2,pos=0)
mtext("H3K4me3",side=2,line=0.5)
mtext("(amplicons I-V)",side=2,line=0.1)

arrows(x[,"days"], x[,"K4.I.ACT2.mean"]-x[,"K4.I.ACT2.SD"], 
       x[,"days"], x[,"K4.I.ACT2.mean"]+x[,"K4.I.ACT2.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.7,1))
arrows(x[,"days"], x[,"K4.II.ACT2.mean"]-x[,"K4.II.ACT2.SD"], 
       x[,"days"], x[,"K4.II.ACT2.mean"]+x[,"K4.II.ACT2.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.2,0.5))
arrows(x[,"days"], x[,"K4.III.ACT2.mean"]-x[,"K4.III.ACT2.SD"], 
       x[,"days"], x[,"K4.III.ACT2.mean"]+x[,"K4.III.ACT2.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,0.2))
arrows(x[,"days"], x[,"K4.IV.ACT2.mean"]-x[,"K4.IV.ACT2.SD"], 
       x[,"days"], x[,"K4.IV.ACT2.mean"]+x[,"K4.IV.ACT2.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col="orange")
arrows(x[,"days"], x[,"K4.V.ACT2.mean"]-x[,"K4.V.ACT2.SD"], 
       x[,"days"], x[,"K4.V.ACT2.mean"]+x[,"K4.V.ACT2.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0.6,0.3,0.6))

arrows(0,1.5,760,1.5,length=0)


######## K4me3(/ACT2) ##############
plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

par(new=T)
plot(x$days,x$K4.VI.ACT2.mean,type="o",pch=16,cex=0.5,lwd=1,
     col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,760),ylim=c(0,0.9))
par(new=T)
plot(x$days,x$K4.VII.ACT2.mean,type="o",pch=18,cex=0.5,lwd=1,
     col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,760),ylim=c(0,0.9))
par(new=T)
plot(x$days,x$K4.VIII.ACT2.mean,type="o",pch=15,cex=0.5,lwd=1,
     col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,760),ylim=c(0,0.9))

axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)

axis(side=2,at=seq(0,0.9,0.3),tcl=0.2,las=2,pos=0)
#mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.5)
#mtext("H3K4me3",side=2,line=0.1)
mtext("H3K4me3",side=2,line=0.5)
mtext("(amplicons VI-VIII)",side=2,line=0.1)

arrows(x$days, x$K4.VI.ACT2.mean-x$K4.VI.ACT2.SD, 
       x$days, x$K4.VI.ACT2.mean+x$K4.VI.ACT2.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.6,0))
arrows(x$days, x$K4.VII.ACT2.mean-x$K4.VII.ACT2.SD, 
       x$days, x$K4.VII.ACT2.mean+x$K4.VII.ACT2.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,1))
arrows(x$days, x$K4.VIII.ACT2.mean-x$K4.VIII.ACT2.SD, 
       x$days, x$K4.VIII.ACT2.mean+x$K4.VIII.ACT2.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0.55,0,0))

arrows(0,0.9,760,0.9,length=0)


######## K27me3(/STM) ##############
plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

x<-read.csv("data/H3K27me3_FLC.csv",header=T,sep=",")

par(new=T)
plot(x[,"days"],x[,"K27.I.STM.mean"],type="o",pch=16,cex=0.5,lwd=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5))
par(new=T)
plot(x[,"days"],x[,"K27.II.STM.mean"],type="o",pch=18,cex=0.5,lwd=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5), lty=1)
par(new=T)
plot(x[,"days"],x[,"K27.III.STM.mean"],type="o",pch=15,cex=0.5,lwd=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5), lty=1)
par(new=T)
plot(x[,"days"],x[,"K27.IV.STM.mean"],type="o",pch=17,cex=0.5,lwd=1,col="orange",axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5), lty=1)
par(new=T)
plot(x[,"days"],x[,"K27.V.STM.mean"],type="o",pch=4,cex=0.5,lwd=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.5), lty=1)

axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1.5,0.5),tcl=0.2,las=2,pos=0)
mtext("H3K27me3",side=2,line=0.5)
mtext("(amplicons I-V)",side=2,line=0.1)

arrows(x[,"days"], x[,"K27.I.STM.mean"]-x[,"K27.I.STM.SD"], 
       x[,"days"], x[,"K27.I.STM.mean"]+x[,"K27.I.STM.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.7,1))
arrows(x[,"days"], x[,"K27.II.STM.mean"]-x[,"K27.II.STM.SD"], 
       x[,"days"], x[,"K27.II.STM.mean"]+x[,"K27.II.STM.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.2,0.5))
arrows(x[,"days"], x[,"K27.III.STM.mean"]-x[,"K27.III.STM.SD"], 
       x[,"days"], x[,"K27.III.STM.mean"]+x[,"K27.III.STM.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,0.2))
arrows(x[,"days"], x[,"K27.IV.STM.mean"]-x[,"K27.IV.STM.SD"], 
       x[,"days"], x[,"K27.IV.STM.mean"]+x[,"K27.IV.STM.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col="orange")
arrows(x[,"days"], x[,"K27.V.STM.mean"]-x[,"K27.V.STM.SD"], 
       x[,"days"], x[,"K27.V.STM.mean"]+x[,"K27.V.STM.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0.6,0.3,0.6))

arrows(0,1.5,760,1.5,length=0)


######## K27me3(/STM) ##############
plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

par(new=T)
plot(x$days,x$K27.VI.STM.mean,type="o",pch=16,cex=0.5,lwd=1,
     col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1))
par(new=T)
plot(x$days,x$K27.VII.STM.mean,type="o",pch=18,cex=0.5,lwd=1,
     col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1), lty=1)
par(new=T)
plot(x$days,x$K27.VIII.STM.mean,type="o",pch=15,cex=0.5,lwd=1,
     col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1), lty=1)

axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1,0.5),tcl=0.2,las=2,pos=0)
#mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.5)
#mtext("H3K27me3",side=2,line=0.1)
mtext("H3K27me3",side=2,line=0.5)
mtext("(amplicons VI-VIII)",side=2,line=0.1)

arrows(x$days, x$K27.VI.STM.mean-x$K27.VI.STM.SD, 
       x$days, x$K27.VI.STM.mean+x$K27.VI.STM.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.6,0))
arrows(x$days, x$K27.VII.STM.mean-x$K27.VII.STM.SD, 
       x$days, x$K27.VII.STM.mean+x$K27.VII.STM.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,1))
arrows(x$days, x$K27.VIII.STM.mean-x$K27.VIII.STM.SD, 
       x$days, x$K27.VIII.STM.mean+x$K27.VIII.STM.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0.55,0,0))

arrows(0,1,760,1,length=0)

set.panel()
dev.off()








pdf("figs/FLC_K27_K4_RNA_FUS3_PP2AA3.pdf",width=3.7,height=4.33)
set.panel(5,1,relax=T)
par(oma=c(0.8,0,0.6,0))
par(mar=c(0.55,2.1,0,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)


###### RNA FLC/PP2AA3 #######################
z<-read.csv("data/120901-140930temp.csv",header=T,sep=",")
z2<-as.matrix(z)

plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

par(new=T)
x<-read.csv("data/RNA_FLC.csv",header=T,sep=",")
x2<-as.matrix(x)
plot(x2[,"days"],x2[,"FLC.PP2AA3.log10.mean"],type="o",pch=16,cex=0.5,
     lwd=1,col=rgb(0,0,0),axes=F,ann=F,xlim=c(0,760),ylim=c(-4,2))
arrows(x2[,"days"], x2[,"FLC.PP2AA3.log10.mean"]-x2[,"FLC.PP2AA3.log10.SD"], 
       x2[,"days"], x2[,"FLC.PP2AA3.log10.mean"]+x2[,"FLC.PP2AA3.log10.SD"], 
       length=.01, angle=90, code=3,lwd=0.5)

axis(side=1,at=place,label=F,tcl=-0.2,pos=-4)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(-4,2,2),lab=c("0.0001","0.01","1","100"),tcl=0.2,las=2,pos=0)
mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.65)
mtext("mRNA",side=2,line=0.3)

arrows(0,2,760,2,length=0)


######## K4me3(/PP2AA3) ##############
plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

x<-read.csv("data/H3K4me3_FLC.csv",header=T,sep=",")

par(new=T)
plot(x[,"days"],x[,"K4.I.PP2AA3.mean"],type="o",pch=16,cex=0.5,lwd=1,
     col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.8))
par(new=T)
plot(x[,"days"],x[,"K4.II.PP2AA3.mean"],type="o",pch=18,cex=0.5,lwd=1,
     col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.8))
par(new=T)
plot(x[,"days"],x[,"K4.III.PP2AA3.mean"],type="o",pch=15,cex=0.5,lwd=1,
     col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.8))
par(new=T)
plot(x[,"days"],x[,"K4.IV.PP2AA3.mean"],type="o",pch=17,cex=0.5,lwd=1,
     col="orange",axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.8))
par(new=T)
plot(x[,"days"],x[,"K4.V.PP2AA3.mean"],type="o",pch=4,cex=0.5,lwd=1,
     col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.8))

axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1.8,0.6),tcl=0.2,las=2,pos=0)
mtext("H3K4me3",side=2,line=0.5)
mtext("(amplicons I-V)",side=2,line=0.1)

arrows(x[,"days"], x[,"K4.I.PP2AA3.mean"]-x[,"K4.I.PP2AA3.SD"], 
       x[,"days"], x[,"K4.I.PP2AA3.mean"]+x[,"K4.I.PP2AA3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.7,1))
arrows(x[,"days"], x[,"K4.II.PP2AA3.mean"]-x[,"K4.II.PP2AA3.SD"], 
       x[,"days"], x[,"K4.II.PP2AA3.mean"]+x[,"K4.II.PP2AA3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.2,0.5))
arrows(x[,"days"], x[,"K4.III.PP2AA3.mean"]-x[,"K4.III.PP2AA3.SD"], 
       x[,"days"], x[,"K4.III.PP2AA3.mean"]+x[,"K4.III.PP2AA3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,0.2))
arrows(x[,"days"], x[,"K4.IV.PP2AA3.mean"]-x[,"K4.IV.PP2AA3.SD"], 
       x[,"days"], x[,"K4.IV.PP2AA3.mean"]+x[,"K4.IV.PP2AA3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col="orange")
arrows(x[,"days"], x[,"K4.V.PP2AA3.mean"]-x[,"K4.V.PP2AA3.SD"], 
       x[,"days"], x[,"K4.V.PP2AA3.mean"]+x[,"K4.V.PP2AA3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0.6,0.3,0.6))

arrows(0,1.8,760,1.8,length=0)


######## K4me3(/PP2AA3) ##############
plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

par(new=T)
plot(x$days,x$K4.VI.PP2AA3.mean,type="o",pch=16,cex=0.5,lwd=1,
     col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.2))
par(new=T)
plot(x$days,x$K4.VII.PP2AA3.mean,type="o",pch=18,cex=0.5,lwd=1,
     col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.2))
par(new=T)
plot(x$days,x$K4.VIII.PP2AA3.mean,type="o",pch=15,cex=0.5,lwd=1,
     col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.2))

axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)

axis(side=2,at=seq(0,1.2,0.4),tcl=0.2,las=2,pos=0)
#mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.5)
#mtext("H3K4me3",side=2,line=0.1)
mtext("H3K4me3",side=2,line=0.5)
mtext("(amplicons VI-VIII)",side=2,line=0.1)

arrows(x$days, x$K4.VI.PP2AA3.mean-x$K4.VI.PP2AA3.SD, 
       x$days, x$K4.VI.PP2AA3.mean+x$K4.VI.PP2AA3.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.6,0))
arrows(x$days, x$K4.VII.PP2AA3.mean-x$K4.VII.PP2AA3.SD, 
       x$days, x$K4.VII.PP2AA3.mean+x$K4.VII.PP2AA3.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,1))
arrows(x$days, x$K4.VIII.PP2AA3.mean-x$K4.VIII.PP2AA3.SD, 
       x$days, x$K4.VIII.PP2AA3.mean+x$K4.VIII.PP2AA3.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0.55,0,0))

arrows(0,1.2,760,1.2,length=0)


######## K27me3(/FUS3) ##############
plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

x<-read.csv("data/H3K27me3_FLC.csv",header=T,sep=",")

par(new=T)
plot(x[,"days"],x[,"K27.I.FUS3.mean"],type="o",pch=16,cex=0.5,lwd=1,
     col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,760),ylim=c(0,2))
par(new=T)
plot(x[,"days"],x[,"K27.II.FUS3.mean"],type="o",pch=18,cex=0.5,lwd=1,
     col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,760),ylim=c(0,2), lty=1)
par(new=T)
plot(x[,"days"],x[,"K27.III.FUS3.mean"],type="o",pch=15,cex=0.5,lwd=1,
     col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,760),ylim=c(0,2), lty=1)
par(new=T)
plot(x[,"days"],x[,"K27.IV.FUS3.mean"],type="o",pch=17,cex=0.5,lwd=1,
     col="orange",axes=F,ann=F,xlim=c(0,760),ylim=c(0,2), lty=1)
par(new=T)
plot(x[,"days"],x[,"K27.V.FUS3.mean"],type="o",pch=4,cex=0.5,lwd=1,
     col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,760),ylim=c(0,2), lty=1)

axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,2,1),tcl=0.2,las=2,pos=0)
mtext("H3K27me3",side=2,line=0.5)
mtext("(amplicons I-V)",side=2,line=0.1)

arrows(x[,"days"], x[,"K27.I.FUS3.mean"]-x[,"K27.I.FUS3.SD"], 
       x[,"days"], x[,"K27.I.FUS3.mean"]+x[,"K27.I.FUS3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.7,1))
arrows(x[,"days"], x[,"K27.II.FUS3.mean"]-x[,"K27.II.FUS3.SD"], 
       x[,"days"], x[,"K27.II.FUS3.mean"]+x[,"K27.II.FUS3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.2,0.5))
arrows(x[,"days"], x[,"K27.III.FUS3.mean"]-x[,"K27.III.FUS3.SD"], 
       x[,"days"], x[,"K27.III.FUS3.mean"]+x[,"K27.III.FUS3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,0.2))
arrows(x[,"days"], x[,"K27.IV.FUS3.mean"]-x[,"K27.IV.FUS3.SD"], 
       x[,"days"], x[,"K27.IV.FUS3.mean"]+x[,"K27.IV.FUS3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col="orange")
arrows(x[,"days"], x[,"K27.V.FUS3.mean"]-x[,"K27.V.FUS3.SD"], 
       x[,"days"], x[,"K27.V.FUS3.mean"]+x[,"K27.V.FUS3.SD"], 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0.6,0.3,0.6))

arrows(0,2,760,2,length=0)


######## K27me3(/FUS3) ##############
plot(z2[,3],type="l",lwd=1, col="gray", axes=F,ann=F,xlim=c(0,760),ylim=c(-5,35))
axis(side=4,at=seq(-5,35,10),tcl=0.2,las=2,pos=760)
mtext("Temperature (°C)",side=4,line=-0.3)

par(new=T)
plot(x$days,x$K27.VI.FUS3.mean,type="o",pch=16,cex=0.5,lwd=1,
     col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.4))
par(new=T)
plot(x$days,x$K27.VII.FUS3.mean,type="o",pch=18,cex=0.5,lwd=1,
     col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.4), lty=1)
par(new=T)
plot(x$days,x$K27.VIII.FUS3.mean,type="o",pch=15,cex=0.5,lwd=1,
     col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,760),ylim=c(0,1.4), lty=1)

axis(side=1,at=place,label=F,tcl=-0.2,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1.4,0.7),tcl=0.2,las=2,pos=0)
#mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.5)
#mtext("H3K27me3",side=2,line=0.1)
mtext("H3K27me3",side=2,line=0.5)
mtext("(amplicons VI-VIII)",side=2,line=0.1)

arrows(x$days, x$K27.VI.FUS3.mean-x$K27.VI.FUS3.SD, 
       x$days, x$K27.VI.FUS3.mean+x$K27.VI.FUS3.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0,0.6,0))
arrows(x$days, x$K27.VII.FUS3.mean-x$K27.VII.FUS3.SD, 
       x$days, x$K27.VII.FUS3.mean+x$K27.VII.FUS3.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(1,0,1))
arrows(x$days, x$K27.VIII.FUS3.mean-x$K27.VIII.FUS3.SD, 
       x$days, x$K27.VIII.FUS3.mean+x$K27.VIII.FUS3.SD, 
       length=.01, angle=90, code=3,lwd=0.5, col=rgb(0.55,0,0))

arrows(0,1.4,760,1.4,length=0)

set.panel()
dev.off()

