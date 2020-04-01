
position <- c(84.5, 313, 1274.5, 2763.5, 4658.5, 5664.5, 6038.5, 6420.5)

######## K27me3(/STM) ##############

pdf("figs/H3K27me3_position.pdf",height=5,width=6)
par(ps=6)
par(mfrow=c(10,5))
par(oma=c(1.4,2,0.6,1))
par(mar=c(0.5,0.7,0.2,0.4))
par(cex=1)
par(mgp=c(0,0.1,0))

x<-read.csv("data/H3K27me3_FLC.csv",header=T,sep=",")
x2<-as.matrix(x)
y<-read.csv("data/H3K27me3_H3K4me3_VI_VII_VIII.csv",header=T,sep=",")
y2<-as.matrix(y)


for(i in 1:45){
  plot(position,c(x2[i,15:19],y2[i,6:8]),type="o",pch=16,cex=0.5,
       axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.5))
  axis(side=1,at=c(0,2000,4000,6000),label=F,tcl=-0.1,pos=0)
  axis(side=2,at=seq(0,1.5,1.5),tcl=0.1,las=2)
  mtext(x2[i,1],side=3,line=-0.6)
  arrows(-500,0,6500,0,length=0)
  #arrows(5500,1.5,5500,0,length=0)
  arrows(position, c(as.numeric(x2[i,15:19]),as.numeric(y2[i,6:8]))-c(as.numeric(x2[i,21:25]),as.numeric(y2[i,18:20])), 
       	 position, c(as.numeric(x2[i,15:19]),as.numeric(y2[i,6:8]))+c(as.numeric(x2[i,21:25]),as.numeric(y2[i,18:20])), 
         length=.01, angle=90, code=3,lwd=0.5)
}

for(i in 46:50){
  plot(position,c(x2[i,15:19],y2[i,6:8]),type="o",pch=16,cex=0.5,
       axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.5))
  axis(side=1,at=c(0,2000,4000,6000),label=F,tcl=-0.1,pos=0)
  axis(side=2,at=seq(0,1.5,1.5),tcl=0.1,las=2)
  mtext(c("0","2,000","4,000","6,000"),at=c(0,2000,4000,6000),side=1,line=-0.3)
  mtext(x2[i,1],side=3,line=-0.6)
  arrows(-500,0,6500,0,length=0)
  #arrows(5500,1.5,5500,0,length=0)
  arrows(position, c(as.numeric(x2[i,15:19]),as.numeric(y2[i,6:8]))-c(as.numeric(x2[i,21:25]),as.numeric(y2[i,18:20])), 
       	 position, c(as.numeric(x2[i,15:19]),as.numeric(y2[i,6:8]))+c(as.numeric(x2[i,21:25]),as.numeric(y2[i,18:20])), 
         length=.01, angle=90, code=3,lwd=0.5)
}


mtext("H3K27me3", side=2, line=0.3, outer=T)
mtext("Distance from TSS (bp)", side=1, line=0, outer=T)
dev.off()





######## K4me3(/ACT2) ##############

pdf("figs/H3K4me3_position.pdf",height=5,width=6)
par(ps=6)
par(mfrow=c(10,5))
par(oma=c(1.4,2,0.6,1))
par(mar=c(0.5,0.7,0.2,0.4))
par(cex=1)
par(mgp=c(0,0.1,0))

x<-read.csv("data/H3K4me3_FLC.csv",header=T,sep=",")
x2<-as.matrix(x)
y<-read.csv("data/H3K27me3_H3K4me3_VI_VII_VIII.csv",header=T,sep=",")
y2<-as.matrix(y)


for(i in 1:45){
  plot(position,c(x2[i,3:7],y2[i,9:11]),type="o",pch=16,cex=0.5,
       axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.5))
  axis(side=1,at=c(0,2000,4000,6000),label=F,tcl=-0.1,pos=0)
  axis(side=2,at=seq(0,1.5,1.5),tcl=0.1,las=2)
  mtext(x2[i,1],side=3,line=-0.6)
  arrows(-500,0,6500,0,length=0)
  #arrows(5500,1.5,5500,0,length=0)
  arrows(position, c(as.numeric(x2[i,3:7]),as.numeric(y2[i,9:11]))-c(as.numeric(x2[i,9:13]),as.numeric(y2[i,21:23])), 
       	 position, c(as.numeric(x2[i,3:7]),as.numeric(y2[i,9:11]))+c(as.numeric(x2[i,9:13]),as.numeric(y2[i,21:23])), 
         length=.01, angle=90, code=3,lwd=0.5)
}

for(i in 46:50){
  plot(position,c(x2[i,3:7],y2[i,9:11]),type="o",pch=16,cex=0.5,
       axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.5))
  axis(side=1,at=c(0,2000,4000,6000),label=F,tcl=-0.1,pos=0)
  axis(side=2,at=seq(0,1.5,1.5),tcl=0.1,las=2)
  mtext(c("0","2,000","4,000","6,000"),at=c(0,2000,4000,6000),side=1,line=-0.3)
  mtext(x2[i,1],side=3,line=-0.6)
  arrows(-500,0,6500,0,length=0)
  #arrows(5500,1.5,5500,0,length=0)
  arrows(position, c(as.numeric(x2[i,3:7]),as.numeric(y2[i,9:11]))-c(as.numeric(x2[i,9:13]),as.numeric(y2[i,21:23])),
       	 position, c(as.numeric(x2[i,3:7]),as.numeric(y2[i,9:11]))+c(as.numeric(x2[i,9:13]),as.numeric(y2[i,21:23])), 
         length=.01, angle=90, code=3,lwd=0.5)
}


mtext("H3K4me3", side=2, line=0.3, outer=T)
mtext("Distance from TSS (bp)", side=1, line=0, outer=T)
dev.off()






##### merged graph #####

######## K27me3(/STM) ##############
pdf("figs/H3K27me3_position_merged_winter.pdf",height=1.2,width=1.8)
par(ps=6)
par(mar=c(1,1.7,1.9,0.8))
par(cex=1)
par(mgp=c(0,0.1,0))
par(xpd=T)

x<-read.csv("data/H3K27me3_FLC.csv",header=T,sep=",")
x2<-as.matrix(x)
y<-read.csv("data/H3K27me3_H3K4me3_VI_VII_VIII.csv",header=T,sep=",")
y2<-as.matrix(y)

plot(position,c(x2[4,15:19],y2[4,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=1)
par(new=T)
plot(position,c(x2[6,15:19],y2[6,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=3)
par(new=T)
plot(position,c(x2[8,15:19],y2[8,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=5)
par(new=T)
plot(position,c(x2[10,15:19],y2[10,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=6)
par(new=T)
plot(position,c(x2[12,15:19],y2[12,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=1,lwd=0.5)
     
for(i in c(4,6,8,10,12)){
	arrows(position, c(as.numeric(x2[i,15:19]),as.numeric(y2[i,6:8]))-c(as.numeric(x2[i,21:25]),as.numeric(y2[i,18:20])), 
    	   position, c(as.numeric(x2[i,15:19]),as.numeric(y2[i,6:8]))+c(as.numeric(x2[i,21:25]),as.numeric(y2[i,18:20])), 
       	   length=.01, angle=90, code=3,lwd=0.5)
}

axis(side=1,at=c(0,2000,4000,6000),label=F,tcl=-0.1,pos=0)
axis(side=2,at=seq(0,1.2,0.6),tcl=0.1,las=2)
mtext(c("0","2,000","4,000","6,000"),at=c(0,2000,4000,6000),side=1,line=-0.45)
#mtext(x2[i,1],side=3,line=-0.6)
arrows(-200,0,6500,0,length=0)
  #arrows(5500,1.5,5500,0,length=0)
mtext("H3K27me3", side=2, line=0.6)
mtext("Distance from TSS (bp)", side=1, line=0.05)

#legend(par()$usr[2]-600, par()$usr[4]+0.12, 
#	c("Nov. 21, 2012","Dec. 18, 2012","Jan. 22, 2013","Feb. 19, 2013","Mar. 19, 2013"),
#	lty=c(1,3,5,6,1), lwd=c(1,1,1,1,0.5),
#	pt.cex=1, y.intersp=0.5, x.intersp=0.1, seg.len=0.7, bty="n")

legend(-2000,2.22,
       legend="Nov. 21, 2012",lty=1,bty="n",x.intersp=0.1,seg.len=1)
legend(-2000,2.02,
       legend="Dec. 18, 2012",lty=3,bty="n",x.intersp=0.1,seg.len=1)
legend(-2000,1.82,
       legend="Jan. 22, 2013",lty=5,bty="n",x.intersp=0.1,seg.len=1)
legend(2100,2.22,
       legend="Feb. 19, 2013",lty=6,bty="n",x.intersp=0.1,seg.len=1)
legend(2100,2.02,
       legend="Mar. 19, 2013",lty=1,lwd=0.5,bty="n",x.intersp=0.1,seg.len=1)
       
dev.off()


pdf("figs/H3K27me3_position_merged_spring.pdf",height=1.2,width=1.8)
par(ps=6)
par(mar=c(1,1.7,1.9,0.8))
par(cex=1)
par(mgp=c(0,0.1,0))
par(xpd=T)

x<-read.csv("data/H3K27me3_FLC.csv",header=T,sep=",")
x2<-as.matrix(x)
y<-read.csv("data/H3K27me3_H3K4me3_VI_VII_VIII.csv",header=T,sep=",")
y2<-as.matrix(y)

plot(position,c(x2[12,15:19],y2[12,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=1)
par(new=T)
plot(position,c(x2[14,15:19],y2[14,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=3)
par(new=T)
plot(position,c(x2[16,15:19],y2[16,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=5)
par(new=T)
plot(position,c(x2[18,15:19],y2[18,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=6)
par(new=T)
plot(position,c(x2[20,15:19],y2[20,6:8]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,1.2),lty=1,lwd=0.5)
     
for(i in c(12,14,16,18,20)){
	arrows(position, c(as.numeric(x2[i,15:19]),as.numeric(y2[i,6:8]))-c(as.numeric(x2[i,21:25]),as.numeric(y2[i,18:20])), 
    	   position, c(as.numeric(x2[i,15:19]),as.numeric(y2[i,6:8]))+c(as.numeric(x2[i,21:25]),as.numeric(y2[i,18:20])), 
       	   length=.01, angle=90, code=3,lwd=0.5)
}

axis(side=1,at=c(0,2000,4000,6000),label=F,tcl=-0.1,pos=0)
axis(side=2,at=seq(0,1.2,0.6),tcl=0.1,las=2)
mtext(c("0","2,000","4,000","6,000"),at=c(0,2000,4000,6000),side=1,line=-0.45)
#mtext(x2[i,1],side=3,line=-0.6)
arrows(-200,0,6500,0,length=0)
#arrows(5500,1.5,5500,0,length=0)
mtext("H3K27me3", side=2, line=0.6)
mtext("Distance from TSS (bp)", side=1, line=0.05)

legend(-2000,2.22,
       legend="Mar. 19, 2013",lty=1,bty="n",x.intersp=0.1,seg.len=1)
legend(-2000,2.02,
       legend="Apr. 16, 2013",lty=3,bty="n",x.intersp=0.1,seg.len=1)
legend(-2000,1.82,
       legend="May 13, 2013",lty=5,bty="n",x.intersp=0.1,seg.len=1)
legend(2100,2.22,
       legend="Jun. 11, 2013",lty=6,bty="n",x.intersp=0.1,seg.len=1)
legend(2100,2.02,
       legend="Jul. 16, 2013",lty=1,lwd=0.5,bty="n",x.intersp=0.1,seg.len=1)

dev.off()





######## K4me3(/ACT2) ##############
pdf("figs/H3K4me3_position_merged_winter.pdf",height=1.2,width=1.8)
par(ps=6)
par(mar=c(1,1.7,1.9,0.8))
par(cex=1)
par(mgp=c(0,0.1,0))
par(xpd=T)

x<-read.csv("data/H3K4me3_FLC.csv",header=T,sep=",")
x2<-as.matrix(x)
y<-read.csv("data/H3K27me3_H3K4me3_VI_VII_VIII.csv",header=T,sep=",")
y2<-as.matrix(y)

plot(position,c(x2[4,3:7],y2[4,9:11]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=1)
par(new=T)
plot(position,c(x2[6,3:7],y2[6,9:11]),type="o",pch=16,cex=0.5,
    axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=3)
par(new=T)
plot(position,c(x2[8,3:7],y2[8,9:11]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=5)
#par(new=T)
#plot(position,c(x2[10,3:7],y2[10,9:11]),type="o",pch=16,cex=0.5,
#     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=6)
par(new=T)
plot(position,c(x2[12,3:7],y2[12,9:11]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=1,lwd=0.5)
     
for(i in c(4,6,8,12)){
	arrows(position, c(as.numeric(x2[i,3:7]),as.numeric(y2[i,9:11]))-c(as.numeric(x2[i,9:13]),as.numeric(y2[i,21:23])), 
    	   position, c(as.numeric(x2[i,3:7]),as.numeric(y2[i,9:11]))+c(as.numeric(x2[i,9:13]),as.numeric(y2[i,21:23])), 
       	   length=.01, angle=90, code=3,lwd=0.5)
}

axis(side=1,at=c(0,2000,4000,6000),label=F,tcl=-0.1,pos=0)
axis(side=2,at=seq(0,0.9,0.3),tcl=0.1,las=2)
mtext(c("0","2,000","4,000","6,000"),at=c(0,2000,4000,6000),side=1,line=-0.45)
#mtext(x2[i,1],side=3,line=-0.6)
arrows(-200,0,6500,0,length=0)
#arrows(5500,1.5,5500,0,length=0)
mtext("H3K4me3", side=2, line=0.6)
mtext("Distance from TSS (bp)", side=1, line=0.05)

legend(-2000,1.665,
       legend="Nov. 21, 2012",lty=1,bty="n",x.intersp=0.1,seg.len=1)
legend(-2000,1.515,
       legend="Dec. 18, 2012",lty=3,bty="n",x.intersp=0.1,seg.len=1)
legend(-2000,1.365,
       legend="Jan. 22, 2013",lty=5,bty="n",x.intersp=0.1,seg.len=1)
legend(2100,1.665,
       legend="Mar. 19, 2013",lty=1,lwd=0.5,bty="n",x.intersp=0.1,seg.len=1)
       
dev.off()


pdf("figs/H3K4me3_position_merged_spring.pdf",height=1.2,width=1.8)
par(ps=6)
par(mar=c(1,1.7,1.9,0.8))
par(cex=1)
par(mgp=c(0,0.1,0))
par(xpd=T)

x<-read.csv("data/H3K4me3_FLC.csv",header=T,sep=",")
x2<-as.matrix(x)
y<-read.csv("data/H3K27me3_H3K4me3_VI_VII_VIII.csv",header=T,sep=",")
y2<-as.matrix(y)

plot(position,c(x2[12,3:7],y2[12,9:11]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=1)
par(new=T)
plot(position,c(x2[14,3:7],y2[14,9:11]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=3)
par(new=T)
plot(position,c(x2[16,3:7],y2[16,9:11]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=5)
#par(new=T)
#plot(position,c(x2[18,3:7],y2[18,9:11]),type="o",pch=16,cex=0.5,
#     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=6)
par(new=T)
plot(position,c(x2[20,3:7],y2[20,9:11]),type="o",pch=16,cex=0.5,
     axes=F,ann=F,xlim=c(0,6500),ylim=c(0,0.9),lty=1,lwd=0.5)
     
for(i in c(12,14,16,20)){
	arrows(position, c(as.numeric(x2[i,3:7]),as.numeric(y2[i,9:11]))-c(as.numeric(x2[i,9:13]),as.numeric(y2[i,21:23])), 
    	   position, c(as.numeric(x2[i,3:7]),as.numeric(y2[i,9:11]))+c(as.numeric(x2[i,9:13]),as.numeric(y2[i,21:23])), 
       	   length=.01, angle=90, code=3,lwd=0.5)
}

axis(side=1,at=c(0,2000,4000,6000),label=F,tcl=-0.1,pos=0)
axis(side=2,at=seq(0,0.9,0.3),tcl=0.1,las=2)
mtext(c("0","2,000","4,000","6,000"),at=c(0,2000,4000,6000),side=1,line=-0.45)
#mtext(x2[i,1],side=3,line=-0.6)
arrows(-200,0,6500,0,length=0)
#arrows(5500,1.5,5500,0,length=0)
mtext("H3K4me3", side=2, line=0.6)
mtext("Distance from TSS (bp)", side=1, line=0.05)

legend(-2000,1.665,
       legend="Mar. 19, 2013",lty=1,bty="n",x.intersp=0.1,seg.len=1)
legend(-2000,1.515,
       legend="Apr. 16, 2013",lty=3,bty="n",x.intersp=0.1,seg.len=1)
legend(-2000,1.365,
       legend="May 13, 2013",lty=5,bty="n",x.intersp=0.1,seg.len=1)
legend(2100,1.665,
       legend="Jul. 16, 2013",lty=1,lwd=0.5,bty="n",x.intersp=0.1,seg.len=1)
       
dev.off()


