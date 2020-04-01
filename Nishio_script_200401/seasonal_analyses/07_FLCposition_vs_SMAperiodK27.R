
######### H3K27me3 ##########################################################
y <- read.csv("data/Histone_rep_SMAmean_120925-140916_log10sma_log10his_result.csv", header=T, sep=",")

period<-c(y$sma[which.max(y$r.sqr.K27.I.STM)],y$sma[which.max(y$r.sqr.K27.II.STM)],
		  y$sma[which.max(y$r.sqr.K27.III.STM)],y$sma[which.max(y$r.sqr.K27.IV.STM)],
		  y$sma[which.max(y$r.sqr.K27.V.STM)],y$sma[which.max(y$r.sqr.K27.VI.STM)],
		  y$sma[which.max(y$r.sqr.K27.VII.STM)],y$sma[which.max(y$r.sqr.K27.VIII.STM)])
position <- c(mean(c(45,124)),mean(c(265,361)),mean(c(1231,1318)),mean(c(2689,2838)),mean(c(4600,4717)),
			  mean(c(5613,5716)),mean(c(5993,6084)),mean(c(6363,6478)))


pdf("figs/SMAperiod-K27_STM.pdf", height=0.9, width=2.7)
par(mar=c(1.2,2.4,1.2,0.2))
par(mgp=c(0,0.15,0))
par(ps=6)
plot(position,period,xlim=c(0,6500),ylim=c(50,150),ylab="",xlab="", yaxt="n",xaxt="n",cex=0.5)
axis(side=1,at=c(0,1000,2000,3000,4000,5000,6000),tcl=-0.1,labels=F)
mtext(c("0","1,000","2,000","3,000","4,000","5,000","6,000"),side=1,
	at=c(0,1000,2000,3000,4000,5000,6000),line=-0.3)
axis(side=2,at=seq(50,150,50),tcl=-0.1,las=2)
mtext(side=1,"Distance from TSS (bp)",line=0.2)
mtext("Moving average", side=2, line=1.25)
mtext("period (days)", side=2, line=0.8)
mtext(expression(paste(italic("AhgFLC"), " H3K27me3")), side=3, line=0)

sp <- smooth.spline(position,period,spar=0.4)
x <- seq(-1000,7000,length=10000)
pred <- predict(sp,x)
lines(pred)
dev.off()

