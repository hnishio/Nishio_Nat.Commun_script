
x<-read.csv("data/Omoide_RNA_histone_1st_year.csv",header=T,sep=",")
x2<-as.matrix(x)
library(calibrate)

nor <- function(x){(x-min(x))/(max(x)-min(x))}

pdf("figs/RNA-K27_Lissajous_ACT2_STM_1st_year.pdf",width=3.3,height=3.5)	
par(ps=6)
par(mfrow=c(3,3))
par(oma=c(1,1,0,0))
par(mar=c(1,1.2,0.9,0.4))
par(cex=1)

### mRNA-H3K27me3 ##################################################

# amplicon I
plot(nor(x2[,"FLC.log10mean"]), nor(x2[,"K27_I_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0,0.7,1),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,24),"FLC.log10mean"]),
       nor(x2[c(1:18,24),"K27_I_STM_mean"]),
       c(1:18,24),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("Amplicon I",side=3,line=0)
abline(1,-1, lwd=1, col="gray")


#amplicon II
plot(nor(x2[,"FLC.log10mean"]), nor(x2[,"K27_II_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0,0.2,0.5),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,24),"FLC.log10mean"]),
       nor(x2[c(1:18,24),"K27_II_STM_mean"]),
       c(1:18,24),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("Amplicon II",side=3,line=0)
abline(1,-1, lwd=1, col="gray")


# amplicon III
plot(nor(x2[,"FLC.log10mean"]), nor(x2[,"K27_III_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(1,0,0.2),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,24),"FLC.log10mean"]),
       nor(x2[c(1:18,24),"K27_III_STM_mean"]),
       c(1:18,24),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("Amplicon III",side=3,line=0)
abline(1,-1, lwd=1, col="gray")


# amplicon IV	
plot(nor(x2[,"FLC.log10mean"]), nor(x2[,"K27_IV_STM_mean"]),type="o",
     pch=16,cex=0.5,col="orange",las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,24),"FLC.log10mean"]),
       nor(x2[c(1:18,24),"K27_IV_STM_mean"]),
       c(1:18,24),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("Amplicon IV",side=3,line=0)
abline(1,-1, lwd=1, col="gray")


# amplicon V
plot(nor(x2[,"FLC.log10mean"]), nor(x2[,"K27_V_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0.6,0.3,0.6),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,24),"FLC.log10mean"]),
       nor(x2[c(1:18,24),"K27_V_STM_mean"]),
       c(1:18,24),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("Amplicon V",side=3,line=0)
abline(1,-1, lwd=1, col="gray")


# amplicon VI
plot(nor(x2[,"FLC.log10mean"]), nor(x2[,"K27_VI_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0,0.6,0),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,24),"FLC.log10mean"]),
       nor(x2[c(1:18,24),"K27_VI_STM_mean"]),
       c(1:18,24),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("Amplicon VI",side=3,line=0)
abline(1,-1, lwd=1, col="gray")


# amplicon VII
plot(nor(x2[,"FLC.log10mean"]), nor(x2[,"K27_VII_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(1,0,1),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,24),"FLC.log10mean"]),
       nor(x2[c(1:18,24),"K27_VII_STM_mean"]),
       c(1:18,24),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("Amplicon VII",side=3,line=0)
abline(1,-1, lwd=1, col="gray")


# amplicon VIII
plot(nor(x2[,"FLC.log10mean"]), nor(x2[,"K27_VIII_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0.55,0,0),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,24),"FLC.log10mean"]),
       nor(x2[c(1:18,24),"K27_VIII_STM_mean"]),
       c(1:18,24),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("Amplicon VIII",side=3,line=0)
abline(1,-1, lwd=1, col="gray")



###mRNA-H3K27me3(All regions)###	

# amplicon I
plot(x2[,"FLC.log10mean"], x2[,"K27_I_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0,0.7,1),las=1,tcl=-0.1, xlim=c(-3.2,0.5),ylim=c(0,1.1),
     xlab="",ylab="",xaxt="n",yaxt="n")

#amplicon II
par(new=T)
plot(x2[,"FLC.log10mean"], x2[,"K27_II_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0,0.2,0.5),las=1,tcl=-0.1, xlim=c(-3.2,0.5),ylim=c(0,1.1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon III
par(new=T)
plot(x2[,"FLC.log10mean"], x2[,"K27_III_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(1,0,0.2),las=1,tcl=-0.1, xlim=c(-3.2,0.5),ylim=c(0,1.1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon IV
par(new=T)
plot(x2[,"FLC.log10mean"], x2[,"K27_IV_STM_mean"],type="o",
     pch=16,cex=0.5,col="orange",las=1,tcl=-0.1, xlim=c(-3.2,0.5),ylim=c(0,1.1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon V
par(new=T)
plot(x2[,"FLC.log10mean"], x2[,"K27_V_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0.6,0.3,0.6),las=1,tcl=-0.1, xlim=c(-3.2,0.5),ylim=c(0,1.1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon VI
par(new=T)
plot(x2[,"FLC.log10mean"], x2[,"K27_VI_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0,0.6,0),las=1,tcl=-0.1, xlim=c(-3.2,0.5),ylim=c(0,1.1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon VII
par(new=T)
plot(x2[,"FLC.log10mean"], x2[,"K27_VII_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(1,0,1),las=1,tcl=-0.1, xlim=c(-3.2,0.5),ylim=c(0,1.1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon VIII
par(new=T)
plot(x2[,"FLC.log10mean"], x2[,"K27_VIII_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0.55,0,0),las=1,tcl=-0.1, xlim=c(-3.2,0.5),ylim=c(0,1.1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)


axis(side=1,at=seq(-3,0,1), tcl=-0.1, labels=F)
mtext(c(expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0)),side=1,at=seq(-3,0,1),line=-0.28)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=0.55)
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=1,line=0.17)
mtext("All regions",side=3,line=0)

dev.off()
