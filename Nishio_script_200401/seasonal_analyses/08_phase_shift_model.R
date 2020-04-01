
curve(sin(x), from=0,to=10)
x<-seq(0,pi*5/2,pi/30)
s<-sin(x)
s1<-sin(x+pi*0/8)
s2<-sin(x+pi*1/8)
s3<-sin(x+pi*2/8)
s4<-sin(x+pi*3/8)
s5<-sin(x+pi*4/8)
s6<-sin(x+pi*5/8)
s7<-sin(x+pi*6/8)
s8<-sin(x+pi*7/8)
s9<-sin(x+pi*8/8)
s10<-sin(x+pi*9/8)
s11<-sin(x+pi*10/8)
s12<-sin(x+pi*11/8)
s13<-sin(x+pi*12/8)
s14<-sin(x+pi*13/8)
s15<-sin(x+pi*14/8)
s16<-sin(x+pi*15/8)
s17<-sin(x+pi*16/8)


pdf("figs/phase_shift_model_190701.pdf",height=4.2, width=4.8)
par(ps=6)
par(mfrow=c(4,5))
par(oma=c(0.4,0.8,0,0))
par(mar=c(0.9,0.7,0.9,0.7))
par(mgp=c(0,0.15,0))
par(cex=1)


### H3K27me3
plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(10,0.3,labels=expression(italic(X)))
text(50-11.25,0.3,labels=expression(italic(Y)),col="red") 
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = -0.75",side=3,line=0)
par(new=T)
plot(s12,type="l",col="red",ylab="",xlab="",axes=F)


plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(10,0.3,labels=expression(italic(X)))
text(50-3.75,0.3,labels=expression(italic(Y)),col="red") 
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = -0.25",side=3,line=0)
par(new=T)
plot(s10,type="l",col="red",ylab="",xlab="",axes=F)


plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(10,0.3,labels=expression(italic(X)))
text(50,0.3,labels=expression(italic(Y)),col="red") 
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = 0.00",side=3,line=0)
par(new=T)
plot(s9,type="l",col="red",ylab="",xlab="",axes=F)


plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(10,0.3,labels=expression(italic(X)))
text(50+3.75,0.3,labels=expression(italic(Y)),col="red") 
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = +0.25",side=3,line=0)
par(new=T)
plot(s8,type="l",col="red",ylab="",xlab="",axes=F)


plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(10,0.3,labels=expression(italic(X)))
text(50+11.25,0.3,labels=expression(italic(Y)),col="red") 
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = +0.75",side=3,line=0)
par(new=T)
plot(s6,type="l",col="red",ylab="",xlab="",axes=F)





plot(s,s12,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)


plot(s,s10,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)


plot(s,s9,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)


plot(s,s8,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)


plot(s,s6,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)



### H3K4me3
plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(27,-0.5,labels=expression(italic(X)))
text(32-11.25,0.5,labels=expression(italic(Y)),col="red")
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = -0.75",side=3,line=0)
par(new=T)
plot(s4,type="l",col="red",ylab="",xlab="",axes=F)


plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(27,-0.5,labels=expression(italic(X)))
text(32-3.75,0.5,labels=expression(italic(Y)),col="red")
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = -0.25",side=3,line=0)
par(new=T)
plot(s2,type="l",col="red",ylab="",xlab="",axes=F)


plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(27,-0.5,labels=expression(italic(X)))
text(32,0.5,labels=expression(italic(Y)),col="red")
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = 0.00",side=3,line=0)
par(new=T)
plot(s17,type="l",col="red",ylab="",xlab="",axes=F)


plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(27,-0.5,labels=expression(italic(X)))
text(32+3.75,0.5,labels=expression(italic(Y)),col="red")
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = +0.25",side=3,line=0)
par(new=T)
plot(s16,type="l",col="red",ylab="",xlab="",axes=F)


plot(s,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(1,76,15),tcl=-0.1,lab=F)
mtext(c("0","1","2","3","4","5"),side=1,at=seq(1,76,15),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,lab=c("0.0","0.5","1.0"),las=1)
text(27,-0.5,labels=expression(italic(X)))
text(32+11.25,0.5,labels=expression(italic(Y)),col="red")
mtext("Level",side=2,line=0.6)
mtext("Time",side=1,line=0.05)
mtext("Phase shift = +0.75",side=3,line=0)
par(new=T)
plot(s14,type="l",col="red",ylab="",xlab="",axes=F)





plot(s,s4,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)


plot(s,s2,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)


plot(s,s17,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)


plot(s,s16,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)


plot(s,s14,type="l",ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=1,at=seq(-1,1,1),tcl=-0.1,lab=F)
mtext(c("0.0","0.5","1.0"),side=1,at=seq(-1,1,1),line=-0.35)
axis(side=2,at=seq(-1,1,1),tcl=-0.1,  pos=-1.08,las=1,lab=c("0.0","0.5","1.0"))
mtext(expression(italic(Y)),side=2,line=0.6)
mtext(expression(italic(X)),side=1,line=0.05)
mtext("V",side=3,line=0,cex=1.5)

dev.off()