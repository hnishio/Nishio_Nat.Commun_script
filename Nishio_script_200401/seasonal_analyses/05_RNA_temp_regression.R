
##### Before analyses, combine RNA data (with replicates) and SMA data #####

RNA <- read.csv("data/RNA_rep.csv",header=T,sep=",")

# daily mean
daily.mean <- read.csv("data/SMAmean_120925-140916.csv",header=T,sep=",")
temp <- NULL
for(i in 1:50){
  temp <- rbind(temp,daily.mean[i,],daily.mean[i,],daily.mean[i,],daily.mean[i,])
}
temp	<- temp+5
RNA.temp <- cbind(RNA,temp)
write.csv(RNA.temp,"data/RNA_rep_SMAmean_120925-140916.csv",quote=F,row.names=F)

# daily max
daily.max <- read.csv("data/SMAmax_120925-140916.csv",header=T,sep=",")
temp <- NULL
for(i in 1:50){
  temp <- rbind(temp,daily.max[i,],daily.max[i,],daily.max[i,],daily.max[i,])
}
temp	<- temp+5
RNA.temp <- cbind(RNA,temp)
write.csv(RNA.temp,"data/RNA_rep_SMAmax_120925-140916.csv",quote=F,row.names=F)

# daily min
daily.min <- read.csv("data/SMAmin_120925-140916.csv",header=T,sep=",")
temp <- NULL
for(i in 1:50){
  temp <- rbind(temp,daily.min[i,],daily.min[i,],daily.min[i,],daily.min[i,])
}
temp	<- temp+5
RNA.temp <- cbind(RNA,temp)
write.csv(RNA.temp,"data/RNA_rep_SMAmin_120925-140916.csv",quote=F,row.names=F)




##### SMA-RNA plot #####
x<-read.csv("data/RNA_rep_SMAmean_120925-140916.csv",header=T,sep=",")
x2<-as.matrix(x)

range(x2[,2])

pdf("figs/Plot_SMAmean-RNA_ACT2.pdf", height=1.1, width=3.5)	
par(ps=6)
par(oma=c(0,0.8,0,0.5))
par(mfrow=c(1,5))
par(mar=c(1.6,0.92,0.8,0))
par(cex=1)

	lapply(list(c(6,2),c(12,2),c(47,2),c(89,2),c(173,2)),
		#1d,1w,6w,12w,24w
		function(z){
     	data<-data.frame(x=x2[,z[1]],y=x2[,z[2]])
	 		model <- lm(y~log10(x), data=data)
	
	 		r2.ca<-round(summary(model)$r.squared,digits=2)
	
			plot(log10(data[,1]),data[,2],las=1,tcl=-0.1,xlim=c(0.7,1.6),ylim=c(-4,2), 
				cex=0.3,ylab="",xlab="", yaxt="n",xaxt="n")
			abline(model)
			text(1.2,1.51,r2.ca)
			text(0.89,1.6,label=expression(paste(R^2,"= ")))
			#axis(side=1,at=seq(0,30,10),tcl=-0.1,las=1,pos=0)

			par(mgp=c(0,0,0))
			axis(side=1,at=seq(1,1.5,0.5),lab=F, tcl=-0.1)
			mtext(c(expression(paste(10^1)),expression(paste(10^1.5))), 
				side=1, at=seq(1,1.5,0.5), line=-0.3)
			mtext("Temperature",side=1,line=0.1)
			mtext("(°C)",side=1,line=0.55)
			
			par(mgp=c(0,0.1,0))
			axis(side=2,at=seq(-4,2,2),lab=F,tcl=-0.1,las=2)
			mtext(c(expression(paste(10^-4)),expression(paste(10^-2)),
					expression(paste(10^0)),expression(paste(10^2))),
				  side=2,line=0.07,at=seq(-4,2,2)+0.2,tcl=-0.1,las=1)
    		    		
    		return(model)
    	}
    )
  
mtext(expression(paste(italic("AhgFLC")," mRNA")),side=2, line=-0.3, outer=T, adj=0.63, cex=1)
mtext("1 day", side=3, line=-0.8, outer=T, adj=0.1, cex=1)
mtext("1 week", side=3, line=-0.8, outer=T, adj=0.31, cex=1)
mtext("6 weeks", side=3, line=-0.8, outer=T, adj=0.53, cex=1)
mtext("12 weeks", side=3, line=-0.8, outer=T, adj=0.75, cex=1)
mtext("24 weeks", side=3, line=-0.8, outer=T, adj=0.98, cex=1)

dev.off()




##### Linear regression between RNA and SMA of temperature #####

# daily mean
rna_sma <- read.csv("data/RNA_rep_SMAmean_120925-140916.csv",header=T,sep=",")

rna_sma<-as.matrix(rna_sma)

sum.table <- matrix(nrow=200,ncol=17)
for (i in 1:200){
	sum.table[i,1] <- i
	
	for (j in 1:4){
	model <- lm(rna_sma[,j+1] ~ log10(rna_sma[,i+5]))
	sum.table[i,(4*j)-3+1] <- summary(model)$coefficient[2,"Estimate"]
	sum.table[i,(4*j)-3+2] <- summary(model)$r.squared
	sum.table[i,(4*j)-3+3] <- AIC(model)
	sum.table[i,(4*j)-3+4] <- logLik(model)
	}
	}
sum.frame <- as.data.frame(sum.table, row.names=F)
names(sum.frame) <- c("sma", 
"coef.RNA.log.ACT2","r.sqr.RNA.log.ACT2","aic.RNA.log.ACT2","loglik.RNA.log.ACT2",
"coef.RNA.ACT2","r.sqr.RNA.ACT2","aic.RNA.ACT2","loglik.RNA.ACT2",
"coef.RNA.log.PP2A","r.sqr.RNA.log.PP2A","aic.RNA.log.PP2A","loglik.RNA.log.PP2A",
"coef.RNA.PP2A","r.sqr.RNA.PP2A","aic.RNA.PP2A","loglik.RNA.PP2A"
)

write.csv(sum.frame,"data/RNA_rep_SMAmean_120925-140916_log10sma_rna_result.csv",quote=F, row.names=F)


# daily max
rna_sma <- read.csv("data/RNA_rep_SMAmax_120925-140916.csv",header=T,sep=",")

rna_sma<-as.matrix(rna_sma)

sum.table <- matrix(nrow=200,ncol=17)
for (i in 1:200){
  sum.table[i,1] <- i
  
  for (j in 1:4){
    model <- lm(rna_sma[,j+1] ~ log10(rna_sma[,i+5]))
    sum.table[i,(4*j)-3+1] <- summary(model)$coefficient[2,"Estimate"]
    sum.table[i,(4*j)-3+2] <- summary(model)$r.squared
    sum.table[i,(4*j)-3+3] <- AIC(model)
    sum.table[i,(4*j)-3+4] <- logLik(model)
  }
}
sum.frame <- as.data.frame(sum.table, row.names=F)
names(sum.frame) <- c("sma", 
                      "coef.RNA.log.ACT2","r.sqr.RNA.log.ACT2","aic.RNA.log.ACT2","loglik.RNA.log.ACT2",
                      "coef.RNA.ACT2","r.sqr.RNA.ACT2","aic.RNA.ACT2","loglik.RNA.ACT2",
                      "coef.RNA.log.PP2A","r.sqr.RNA.log.PP2A","aic.RNA.log.PP2A","loglik.RNA.log.PP2A",
                      "coef.RNA.PP2A","r.sqr.RNA.PP2A","aic.RNA.PP2A","loglik.RNA.PP2A"
)

write.csv(sum.frame,"data/RNA_rep_SMAmax_120925-140916_log10sma_rna_result.csv",quote=F, row.names=F)



# daily min
rna_sma <- read.csv("data/RNA_rep_SMAmin_120925-140916.csv",header=T,sep=",")

rna_sma<-as.matrix(rna_sma)

sum.table <- matrix(nrow=200,ncol=17)
for (i in 1:200){
  sum.table[i,1] <- i
  
  for (j in 1:4){
    model <- lm(rna_sma[,j+1] ~ log10(rna_sma[,i+5]))
    sum.table[i,(4*j)-3+1] <- summary(model)$coefficient[2,"Estimate"]
    sum.table[i,(4*j)-3+2] <- summary(model)$r.squared
    sum.table[i,(4*j)-3+3] <- AIC(model)
    sum.table[i,(4*j)-3+4] <- logLik(model)
  }
}
sum.frame <- as.data.frame(sum.table, row.names=F)
names(sum.frame) <- c("sma", 
                      "coef.RNA.log.ACT2","r.sqr.RNA.log.ACT2","aic.RNA.log.ACT2","loglik.RNA.log.ACT2",
                      "coef.RNA.ACT2","r.sqr.RNA.ACT2","aic.RNA.ACT2","loglik.RNA.ACT2",
                      "coef.RNA.log.PP2A","r.sqr.RNA.log.PP2A","aic.RNA.log.PP2A","loglik.RNA.log.PP2A",
                      "coef.RNA.PP2A","r.sqr.RNA.PP2A","aic.RNA.PP2A","loglik.RNA.PP2A"
)

write.csv(sum.frame,"data/RNA_rep_SMAmin_120925-140916_log10sma_rna_result.csv",quote=F, row.names=F)




##### Graph of R2 value #####

# daily mean
y <- read.csv("data/RNA_rep_SMAmean_120925-140916_log10sma_rna_result.csv", header=T, sep=",")
y<-as.matrix(y)

# FLC logRNA /ACT2
pdf("figs/R2_SMAmean-logRNA_ACT2.pdf", height=1.3, width=1.38)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,3],type="l",pch=16,cex=1,axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " mRNA")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)
dev.off()

y[,1][which.max(y[,3])]


# FLC logRNA /PP2AA3
pdf("figs/R2_SMAmean-logRNA_PP2AA3.pdf", height=1.3, width=1.38)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,11],type="l",pch=16,cex=1,axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " mRNA (/", italic("AhgPP2AA3"), ")")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)
dev.off()

y[,1][which.max(y[,11])]




# daily max
y <- read.csv("data/RNA_rep_SMAmax_120925-140916_log10sma_rna_result.csv", header=T, sep=",")
y<-as.matrix(y)

# FLC logRNA /ACT2
pdf("figs/R2_SMAmax-logRNA_ACT2.pdf", height=1.3, width=1.38)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,3],type="l",pch=16,cex=1,axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " mRNA")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)
dev.off()

y[,1][which.max(y[,3])]





# daily min
y <- read.csv("data/RNA_rep_SMAmin_120925-140916_log10sma_rna_result.csv", header=T, sep=",")
y<-as.matrix(y)

# FLC logRNA /ACT2
pdf("figs/R2_SMAmin-logRNA_ACT2.pdf", height=1.3, width=1.38)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,3],type="l",pch=16,cex=1,axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " mRNA")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)
dev.off()

y[,1][which.max(y[,3])]


