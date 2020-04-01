
x<-read.csv("data/Omoide_RNA_histone_2nd_year.csv",header=T,sep=",")
x2<-as.matrix(x)
library(calibrate)

nor <- function(x){(x-min(x))/(max(x)-min(x))}

pdf("figs/K4-K27_Lissajous_ACT2_STM_2nd_year.pdf",width=3.4,height=3.5)	
par(ps=6)
par(mfrow=c(3,3))
par(oma=c(1,0,0,0))
par(mar=c(1,1.7,0.9,0.4))
par(cex=1)

### H3K4me3-H3K27me3 ##################################################

# amplicon I
plot(nor(x2[,"K4_I_ACT2_mean"]), nor(x2[,"K27_I_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0,0.7,1),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,26),"K4_I_ACT2_mean"]),
       nor(x2[c(1:18,26),"K27_I_STM_mean"]),
       c(1:18,26),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(amplicon I)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)
abline(1,-1, lwd=1, col="gray")


#amplicon II
plot(nor(x2[,"K4_I_ACT2_mean"]), nor(x2[,"K27_II_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0,0.2,0.5),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,26),"K4_I_ACT2_mean"]),
       nor(x2[c(1:18,26),"K27_II_STM_mean"]),
       c(1:18,26),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(amplicon II)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)
abline(1,-1, lwd=1, col="gray")


# amplicon III
plot(nor(x2[,"K4_I_ACT2_mean"]), nor(x2[,"K27_III_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(1,0,0.2),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,26),"K4_I_ACT2_mean"]),
       nor(x2[c(1:18,26),"K27_III_STM_mean"]),
       c(1:18,26),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(amplicon III)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)
abline(1,-1, lwd=1, col="gray")


# amplicon IV	
plot(nor(x2[,"K4_I_ACT2_mean"]), nor(x2[,"K27_IV_STM_mean"]),type="o",
     pch=16,cex=0.5,col="orange",las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,26),"K4_I_ACT2_mean"]),
       nor(x2[c(1:18,26),"K27_IV_STM_mean"]),
       c(1:18,26),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(amplicon IV)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)
abline(1,-1, lwd=1, col="gray")


# amplicon V
plot(nor(x2[,"K4_I_ACT2_mean"]), nor(x2[,"K27_V_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0.6,0.3,0.6),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,26),"K4_I_ACT2_mean"]),
       nor(x2[c(1:18,26),"K27_V_STM_mean"]),
       c(1:18,26),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(amplicon V)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)
abline(1,-1, lwd=1, col="gray")


# amplicon VI
plot(nor(x2[,"K4_I_ACT2_mean"]), nor(x2[,"K27_VI_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0,0.6,0),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,26),"K4_I_ACT2_mean"]),
       nor(x2[c(1:18,26),"K27_VI_STM_mean"]),
       c(1:18,26),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(amplicon VI)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)
abline(1,-1, lwd=1, col="gray")


# amplicon VII
plot(nor(x2[,"K4_I_ACT2_mean"]), nor(x2[,"K27_VII_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(1,0,1),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,26),"K4_I_ACT2_mean"]),
       nor(x2[c(1:18,26),"K27_VII_STM_mean"]),
       c(1:18,26),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(amplicon VII)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)
abline(1,-1, lwd=1, col="gray")


# amplicon VIII
plot(nor(x2[,"K4_I_ACT2_mean"]), nor(x2[,"K27_VIII_STM_mean"]),type="o",
     pch=16,cex=0.5,col=rgb(0.55,0,0),las=1,tcl=-0.1, xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

par(ps=5)
textxy(nor(x2[c(1:18,26),"K4_I_ACT2_mean"]),
       nor(x2[c(1:18,26),"K27_VIII_STM_mean"]),
       c(1:18,26),cex=1)	

par(ps=6)
par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(amplicon VIII)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)
abline(1,-1, lwd=1, col="gray")



###H3K4me3-H3K27me3(All regions)###	

# amplicon I
plot(x2[,"K4_I_ACT2_mean"], x2[,"K27_I_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0,0.7,1),las=1,tcl=-0.1, xlim=c(0,1.1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n")

#amplicon II
par(new=T)
plot(x2[,"K4_I_ACT2_mean"], x2[,"K27_II_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0,0.2,0.5),las=1,tcl=-0.1, xlim=c(0,1.1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon III
par(new=T)
plot(x2[,"K4_I_ACT2_mean"], x2[,"K27_III_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(1,0,0.2),las=1,tcl=-0.1, xlim=c(0,1.1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon IV
par(new=T)
plot(x2[,"K4_I_ACT2_mean"], x2[,"K27_IV_STM_mean"],type="o",
     pch=16,cex=0.5,col="orange",las=1,tcl=-0.1, xlim=c(0,1.1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon V
par(new=T)
plot(x2[,"K4_I_ACT2_mean"], x2[,"K27_V_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0.6,0.3,0.6),las=1,tcl=-0.1, xlim=c(0,1.1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon VI
par(new=T)
plot(x2[,"K4_I_ACT2_mean"], x2[,"K27_VI_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0,0.6,0),las=1,tcl=-0.1, xlim=c(0,1.1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon VII
par(new=T)
plot(x2[,"K4_I_ACT2_mean"], x2[,"K27_VII_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(1,0,1),las=1,tcl=-0.1, xlim=c(0,1.1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

# amplicon VIII
par(new=T)
plot(x2[,"K4_I_ACT2_mean"], x2[,"K27_VIII_STM_mean"],type="o",
     pch=16,cex=0.5,col=rgb(0.55,0,0),las=1,tcl=-0.1, xlim=c(0,1.1),ylim=c(0,1),
     xlab="",ylab="",xaxt="n",yaxt="n",axes=F, ann=F)

par(mgp=c(0,-0.35,0))
axis(side=1,at=seq(0,1,0.5), tcl=-0.1)
par(mgp=c(0,0.15,0))
axis(side=2,at=seq(0,1,0.5), tcl=-0.1, las=2)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2,line=1)
mtext("(all regions)",side=2,line=0.65)
mtext(expression(paste(italic("AhgFLC")," H3K4me3")),side=1,line=0.17)
mtext("(amplicon I)",side=1,line=0.52)

dev.off()
