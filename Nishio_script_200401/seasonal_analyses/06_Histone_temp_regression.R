
##### Before analyses, combine histone data (with replicates) and SMA data #####

histone <- read.csv("data/Histone_rep.csv",header=T,sep=",")

# daily mean
daily.mean <- read.csv("data/SMAmean_120925-140916.csv",header=T,sep=",")
temp <- NULL
for(i in 1:50){
  temp <- rbind(temp,daily.mean[i,],daily.mean[i,],daily.mean[i,],daily.mean[i,])
}
temp	<- temp+5
his.temp <- cbind(histone,temp)
write.csv(his.temp,"data/Histone_rep_SMAmean_120925-140916.csv",quote=F,row.names=F)

# daily max
daily.max <- read.csv("data/SMAmax_120925-140916.csv",header=T,sep=",")
temp <- NULL
for(i in 1:50){
  temp <- rbind(temp,daily.max[i,],daily.max[i,],daily.max[i,],daily.max[i,])
}
temp	<- temp+5
his.temp <- cbind(histone,temp)
write.csv(his.temp,"data/Histone_rep_SMAmax_120925-140916.csv",quote=F,row.names=F)

# daily min
daily.min <- read.csv("data/SMAmin_120925-140916.csv",header=T,sep=",")
temp <- NULL
for(i in 1:50){
  temp <- rbind(temp,daily.min[i,],daily.min[i,],daily.min[i,],daily.min[i,])
}
temp	<- temp+5
his.temp <- cbind(histone,temp)
write.csv(his.temp,"data/Histone_rep_SMAmin_120925-140916.csv",quote=F,row.names=F)



# Check of the minimum SMA
daily.mean2 <- c(as.matrix(daily.mean))
daily.mean2[which(daily.mean2<0)]
daily.mean2[order(daily.mean2)][1:10]

daily.max2 <- c(as.matrix(daily.max))
daily.max2[which(daily.max2<0)]
daily.max2[order(daily.max2)][1]

daily.min2 <- c(as.matrix(daily.min))
daily.min2[which(daily.min2<0)]
length(daily.min2[which(daily.min2<0)])
daily.min2[order(daily.min2)][1]



##### SMA-K27 plot #####
x<-read.csv("data/Histone_rep_SMAmean_120925-140916.csv",header=T,sep=",")
x2<-as.matrix(x)

max(log10(x2[,40]))

###### H3K27me3(log10-log10) STM #########
pdf("figs/Plot_SMAmean-K27_STM.pdf", height=1.1, width=3.5)	
par(ps=6)
par(oma=c(0,0.8,0,0.5))
par(mfrow=c(1,5))
par(mar=c(1.6,0.92,0.8,0))
par(cex=1)

	lapply(list(c(34,10),c(40,10),c(75,10),c(117,10),c(201,10)),
		#1d,1w,6w,12w,24w
		function(z){
     	data<-data.frame(x=x2[,z[1]],y=x2[,z[2]])
	 		model <- lm(log10(y)~log10(x), data=data)
	
	 		r2.ca<-round(summary(model)$r.squared,digits=2)
	
			plot(log10(data[,1]),log10(data[,2]),las=1,tcl=-0.1,xlim=c(0.7,1.6),ylim=c(-3,1), 
				cex=0.3,ylab="",xlab="", yaxt="n",xaxt="n",col=rgb(0,0.7,1))
			abline(model)
			text(1.43,0.71,r2.ca)
			text(1.12,0.8,label=expression(paste(R^2,"= ")))
			#axis(side=1,at=seq(0,30,10),tcl=-0.1,las=1,pos=0)

			par(mgp=c(0,0,0))
			axis(side=1,at=seq(1,1.5,0.5),lab=F, tcl=-0.1)
			mtext(c(expression(paste(10^1)),expression(paste(10^1.5))), 
				side=1, at=seq(1,1.5,0.5), line=-0.3)
			mtext("Temperature",side=1,line=0.1)
			mtext("(°C)",side=1,line=0.55)
			
			par(mgp=c(0,0.1,0))
			axis(side=2,at=seq(-3,1,1),lab=F,tcl=-0.1,las=2)
			mtext(c(expression(paste(10^-3)),expression(paste(10^-2)),expression(paste(10^-1)),
				expression(paste(10^0)),expression(paste(10^1))),side=2,line=0.07,at=seq(-2.9,1.1,1),tcl=-0.1,las=1)
    		    		
    		return(model)
    	}
    )
  
mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2, line=-0.3, outer=T, adj=0.63, cex=1)
mtext("1 day", side=3, line=-0.8, outer=T, adj=0.1, cex=1)
mtext("1 week", side=3, line=-0.8, outer=T, adj=0.31, cex=1)
mtext("6 weeks", side=3, line=-0.8, outer=T, adj=0.53, cex=1)
mtext("12 weeks", side=3, line=-0.8, outer=T, adj=0.75, cex=1)
mtext("24 weeks", side=3, line=-0.8, outer=T, adj=0.98, cex=1)

dev.off()




##### SMA-K27 plot #####
x<-read.csv("data/Histone_rep_SMAmin_120925-140916.csv",header=T,sep=",")
x2<-as.matrix(x)

###### H3K27me3(log10-log10) STM #########
pdf("figs/Plot_SMAmin-K27_STM.pdf", height=1.1, width=3.5)	
par(ps=6)
par(oma=c(0,0.8,0,0.5))
par(mfrow=c(1,5))
par(mar=c(1.6,0.92,0.8,0))
par(cex=1)

lapply(list(c(34,10),c(40,10),c(75,10),c(117,10),c(201,10)),
       #1d,1w,6w,12w,24w
       function(z){
         data<-data.frame(x=x2[,z[1]],y=x2[,z[2]])
         model <- lm(log10(y)~log10(x), data=data)
         
         r2.ca<-round(summary(model)$r.squared,digits=2)
         
         plot(log10(data[,1]),log10(data[,2]),las=1,tcl=-0.1,xlim=c(-1.5,1.5),ylim=c(-3,1), 
              cex=0.3,ylab="",xlab="", yaxt="n",xaxt="n",col=rgb(0,0.7,1))
         abline(model)
         text(1.2,0.71,r2.ca)
         text(0.7,0.8,label=expression(paste(R^2,"= ")))
         #axis(side=1,at=seq(0,30,10),tcl=-0.1,las=1,pos=0)
         
         par(mgp=c(0,0,0))
         axis(side=1,at=seq(-1.5,1.5,0.5),lab=F, tcl=-0.1)
         mtext(c(expression(paste(10^0)),"",expression(paste(10^1)),""), 
               side=1, at=c(0,0.5,1,1.5), line=-0.3)
         mtext("Temperature",side=1,line=0.1)
         mtext("(°C)",side=1,line=0.55)
         
         par(mgp=c(0,0.1,0))
         axis(side=2,at=seq(-3,1,1),lab=F,tcl=-0.1,las=2)
         mtext(c(expression(paste(10^-3)),expression(paste(10^-2)),expression(paste(10^-1)),
                 expression(paste(10^0)),expression(paste(10^1))),side=2,line=0.07,at=seq(-2.9,1.1,1),tcl=-0.1,las=1)
         
         return(model)
       }
)

mtext(expression(paste(italic("AhgFLC")," H3K27me3")),side=2, line=-0.3, outer=T, adj=0.63, cex=1)
mtext("1 day", side=3, line=-0.8, outer=T, adj=0.1, cex=1)
mtext("1 week", side=3, line=-0.8, outer=T, adj=0.31, cex=1)
mtext("6 weeks", side=3, line=-0.8, outer=T, adj=0.53, cex=1)
mtext("12 weeks", side=3, line=-0.8, outer=T, adj=0.75, cex=1)
mtext("24 weeks", side=3, line=-0.8, outer=T, adj=0.98, cex=1)

dev.off()





##### Linear regression against SMAs with one-day differences in window length #####

# daily mean
his.sma <- read.csv("data/Histone_rep_SMAmean_120925-140916.csv",header=T,sep=",")
his.sma <- as.matrix(his.sma)

sum.table <- matrix(nrow=200,ncol=129)
for (i in 1:200){
	sum.table[i,1] <- i
	
	for (j in 1:32){
	model <- lm(log10(his.sma[,j+1]) ~ log10(his.sma[,i+33]))
	sum.table[i,(4*j)-3+1] <- summary(model)$coefficient[2,"Estimate"]
	sum.table[i,(4*j)-3+2] <- summary(model)$r.squared
	sum.table[i,(4*j)-3+3] <- AIC(model)
	sum.table[i,(4*j)-3+4] <- logLik(model)
	}
	}
sum.frame <- as.data.frame(sum.table, row.names=F)

names(sum.frame) <-
  c("sma",
    paste0(c("coef.","r.sqr.","aic.","loglik."), 
           c(rep("K27.I.FUS3",4),rep("K27.II.FUS3",4),rep("K27.III.FUS3",4),rep("K27.IV.FUS3",4),
             rep("K27.V.FUS3",4),rep("K27.VI.FUS3",4),rep("K27.VII.FUS3",4),rep("K27.VIII.FUS3",4),
             rep("K27.I.STM",4),rep("K27.II.STM",4),rep("K27.III.STM",4),rep("K27.IV.STM",4),
             rep("K27.V.STM",4),rep("K27.VI.STM",4),rep("K27.VII.STM",4),rep("K27.VIII.STM",4),          
             rep("K4.I.ACT2",4),rep("K4.II.ACT2",4),rep("K4.III.ACT2",4),rep("K4.IV.ACT2",4),
             rep("K4.V.ACT2",4),rep("K4.VI.ACT2",4),rep("K4.VII.ACT2",4),rep("K4.VIII.ACT2",4),             
             rep("K4.I.PP2AA3",4),rep("K4.II.PP2AA3",4),rep("K4.III.PP2AA3",4),rep("K4.IV.PP2AA3",4),
             rep("K4.V.PP2AA3",4),rep("K4.VI.PP2AA3",4),rep("K4.VII.PP2AA3",4),rep("K4.VIII.PP2AA3",4)           
          ))
  )

write.csv(sum.frame,"data/Histone_rep_SMAmean_120925-140916_log10sma_log10his_result.csv",quote=F,row.names=F)




# daily max
his.sma <- read.csv("data/Histone_rep_SMAmax_120925-140916.csv",header=T,sep=",")
his.sma <- as.matrix(his.sma)

sum.table <- matrix(nrow=200,ncol=129)
for (i in 1:200){
  sum.table[i,1] <- i
  
  for (j in 1:32){
    model <- lm(log10(his.sma[,j+1]) ~ log10(his.sma[,i+33]))
    sum.table[i,(4*j)-3+1] <- summary(model)$coefficient[2,"Estimate"]
    sum.table[i,(4*j)-3+2] <- summary(model)$r.squared
    sum.table[i,(4*j)-3+3] <- AIC(model)
    sum.table[i,(4*j)-3+4] <- logLik(model)
  }
}
sum.frame <- as.data.frame(sum.table, row.names=F)

names(sum.frame) <-
  c("sma",
    paste0(c("coef.","r.sqr.","aic.","loglik."), 
           c(rep("K27.I.FUS3",4),rep("K27.II.FUS3",4),rep("K27.III.FUS3",4),rep("K27.IV.FUS3",4),
             rep("K27.V.FUS3",4),rep("K27.VI.FUS3",4),rep("K27.VII.FUS3",4),rep("K27.VIII.FUS3",4),
             rep("K27.I.STM",4),rep("K27.II.STM",4),rep("K27.III.STM",4),rep("K27.IV.STM",4),
             rep("K27.V.STM",4),rep("K27.VI.STM",4),rep("K27.VII.STM",4),rep("K27.VIII.STM",4),          
             rep("K4.I.ACT2",4),rep("K4.II.ACT2",4),rep("K4.III.ACT2",4),rep("K4.IV.ACT2",4),
             rep("K4.V.ACT2",4),rep("K4.VI.ACT2",4),rep("K4.VII.ACT2",4),rep("K4.VIII.ACT2",4),             
             rep("K4.I.PP2AA3",4),rep("K4.II.PP2AA3",4),rep("K4.III.PP2AA3",4),rep("K4.IV.PP2AA3",4),
             rep("K4.V.PP2AA3",4),rep("K4.VI.PP2AA3",4),rep("K4.VII.PP2AA3",4),rep("K4.VIII.PP2AA3",4)           
           ))
  )

write.csv(sum.frame,"data/Histone_rep_SMAmax_120925-140916_log10sma_log10his_result.csv",quote=F,row.names=F)



# daily min
his.sma <- read.csv("data/Histone_rep_SMAmin_120925-140916.csv",header=T,sep=",")
his.sma <- as.matrix(his.sma)

sum.table <- matrix(nrow=200,ncol=129)
for (i in 1:200){
  sum.table[i,1] <- i
  
  for (j in 1:32){
    model <- lm(log10(his.sma[,j+1]) ~ log10(his.sma[,i+33]))
    sum.table[i,(4*j)-3+1] <- summary(model)$coefficient[2,"Estimate"]
    sum.table[i,(4*j)-3+2] <- summary(model)$r.squared
    sum.table[i,(4*j)-3+3] <- AIC(model)
    sum.table[i,(4*j)-3+4] <- logLik(model)
  }
}
sum.frame <- as.data.frame(sum.table, row.names=F)

names(sum.frame) <-
  c("sma",
    paste0(c("coef.","r.sqr.","aic.","loglik."), 
           c(rep("K27.I.FUS3",4),rep("K27.II.FUS3",4),rep("K27.III.FUS3",4),rep("K27.IV.FUS3",4),
             rep("K27.V.FUS3",4),rep("K27.VI.FUS3",4),rep("K27.VII.FUS3",4),rep("K27.VIII.FUS3",4),
             rep("K27.I.STM",4),rep("K27.II.STM",4),rep("K27.III.STM",4),rep("K27.IV.STM",4),
             rep("K27.V.STM",4),rep("K27.VI.STM",4),rep("K27.VII.STM",4),rep("K27.VIII.STM",4),          
             rep("K4.I.ACT2",4),rep("K4.II.ACT2",4),rep("K4.III.ACT2",4),rep("K4.IV.ACT2",4),
             rep("K4.V.ACT2",4),rep("K4.VI.ACT2",4),rep("K4.VII.ACT2",4),rep("K4.VIII.ACT2",4),             
             rep("K4.I.PP2AA3",4),rep("K4.II.PP2AA3",4),rep("K4.III.PP2AA3",4),rep("K4.IV.PP2AA3",4),
             rep("K4.V.PP2AA3",4),rep("K4.VI.PP2AA3",4),rep("K4.VII.PP2AA3",4),rep("K4.VIII.PP2AA3",4)           
           ))
  )

write.csv(sum.frame,"data/Histone_rep_SMAmin_120925-140916_log10sma_log10his_result.csv",quote=F,row.names=F)




##### Graph of R2 value #####

# daily mean
y <- read.csv("data/Histone_rep_SMAmean_120925-140916_log10sma_log10his_result.csv", header=T, sep=",")
y<-as.matrix(y)

# K27_FUS3
pdf("figs/R2_SMAmean-K27_FUS3_all.pdf", height=1.3, width=1.4)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,3],type="l",pch=16,cex=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,7],type="l",pch=18,cex=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,11],type="l",pch=15,cex=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,15],type="l",pch=17,cex=1,col="orange",axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,19],type="l",pch=18,cex=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,23],type="l",pch=18,cex=1,col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,27],type="l",pch=18,cex=1,col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,31],type="l",pch=18,cex=1,col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " H3K27me3 (/", italic("AhgFUS3"), ")")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)

legend(par()$usr[2]-35, par()$usr[4]+0.12, c("I","II","III","IV","V","VI","VII","VIII"), lwd=1,
col=c(rgb(0,0.7,1), rgb(0,0.2,0.5), rgb(1,0,0.2), "orange", rgb(0.6,0.3,0.6), rgb(0,0.6,0), rgb(1,0,1), rgb(0.55,0,0)),lty=c(1,1,1,1,1,6,6,6),
pt.cex=1, y.intersp=0.45, x.intersp=0.1, seg.len=0.7, bty="n")
dev.off()


# K27_STM
pdf("figs/R2_SMAmean-K27_STM_all.pdf", height=1.3, width=1.4)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,35],type="l",pch=16,cex=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,39],type="l",pch=18,cex=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,43],type="l",pch=15,cex=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,47],type="l",pch=17,cex=1,col="orange",axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,51],type="l",pch=18,cex=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,55],type="l",pch=18,cex=1,col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,59],type="l",pch=18,cex=1,col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,63],type="l",pch=18,cex=1,col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " H3K27me3")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)

legend(par()$usr[2]-35, par()$usr[4]+0.12, c("I","II","III","IV","V","VI","VII","VIII"), lwd=1,
col=c(rgb(0,0.7,1), rgb(0,0.2,0.5), rgb(1,0,0.2), "orange", rgb(0.6,0.3,0.6), rgb(0,0.6,0), rgb(1,0,1), rgb(0.55,0,0)),lty=c(1,1,1,1,1,6,6,6),
pt.cex=1, y.intersp=0.45, x.intersp=0.1, seg.len=0.7, bty="n")
dev.off()

y[,1][which.max(y[,35])]
y[,1][which.max(y[,39])]
y[,1][which.max(y[,43])]
y[,1][which.max(y[,47])]
y[,1][which.max(y[,51])]
y[,1][which.max(y[,55])]
y[,1][which.max(y[,59])]
y[,1][which.max(y[,63])]


# K4_ACT2
pdf("figs/R2_SMAmean-K4_ACT2_all.pdf", height=1.3, width=1.4)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,67],type="l",pch=16,cex=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,71],type="l",pch=18,cex=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,75],type="l",pch=15,cex=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,79],type="l",pch=17,cex=1,col="orange",axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,83],type="l",pch=18,cex=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,87],type="l",pch=18,cex=1,col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,91],type="l",pch=18,cex=1,col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,95],type="l",pch=18,cex=1,col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " H3K4me3")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)

legend(90, 1.14, c("I","II","VI","VII","VIII"), lwd=1,
col=c(rgb(0,0.7,1), rgb(0,0.2,0.5), rgb(0,0.6,0), rgb(1,0,1), rgb(0.55,0,0)),lty=c(1,1,6,6,6),
pt.cex=1, y.intersp=0.45, x.intersp=0.1, seg.len=0.7, bty="n")
dev.off()


y[,1][which.max(y[,67])]
y[,1][which.max(y[,71])]
#y[,1][which.max(y[,75])]
#y[,1][which.max(y[,79])]
#y[,1][which.max(y[,83])]
y[,1][which.max(y[,87])]
y[,1][which.max(y[,91])]
y[,1][which.max(y[,95])]


# K4_PP2AA3
pdf("figs/R2_SMAmean-K4_PP2AA3_all.pdf", height=1.3, width=1.4)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,99],type="l",pch=16,cex=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,103],type="l",pch=18,cex=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,107],type="l",pch=15,cex=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,111],type="l",pch=17,cex=1,col="orange",axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,115],type="l",pch=18,cex=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,119],type="l",pch=18,cex=1,col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,123],type="l",pch=18,cex=1,col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,127],type="l",pch=18,cex=1,col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " H3K4me3 (/", italic("AhgPP2AA3"), ")")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)

legend(90, 1.14, c("I","II","VI","VII","VIII"), lwd=1,
col=c(rgb(0,0.7,1), rgb(0,0.2,0.5), rgb(0,0.6,0), rgb(1,0,1), rgb(0.55,0,0)),lty=c(1,1,6,6,6),
pt.cex=1, y.intersp=0.45, x.intersp=0.1, seg.len=0.7, bty="n")

dev.off()






# daily max
y <- read.csv("data/Histone_rep_SMAmax_120925-140916_log10sma_log10his_result.csv", header=T, sep=",")
y<-as.matrix(y)

# K27_STM
pdf("figs/R2_SMAmax-K27_STM_all.pdf", height=1.3, width=1.4)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,35],type="l",pch=16,cex=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,39],type="l",pch=18,cex=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,43],type="l",pch=15,cex=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,47],type="l",pch=17,cex=1,col="orange",axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,51],type="l",pch=18,cex=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,55],type="l",pch=18,cex=1,col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,59],type="l",pch=18,cex=1,col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,63],type="l",pch=18,cex=1,col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " H3K27me3")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)

legend(par()$usr[2]-35, par()$usr[4]+0.12, c("I","II","III","IV","V","VI","VII","VIII"), lwd=1,
       col=c(rgb(0,0.7,1), rgb(0,0.2,0.5), rgb(1,0,0.2), "orange", rgb(0.6,0.3,0.6), rgb(0,0.6,0), rgb(1,0,1), rgb(0.55,0,0)),lty=c(1,1,1,1,1,6,6,6),
       pt.cex=1, y.intersp=0.45, x.intersp=0.1, seg.len=0.7, bty="n")
dev.off()

y[,1][which.max(y[,35])]
y[,1][which.max(y[,39])]
y[,1][which.max(y[,43])]
y[,1][which.max(y[,47])]
y[,1][which.max(y[,51])]
y[,1][which.max(y[,55])]
y[,1][which.max(y[,59])]
y[,1][which.max(y[,63])]


# K4_ACT2
pdf("figs/R2_SMAmax-K4_ACT2_all.pdf", height=1.3, width=1.4)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,67],type="l",pch=16,cex=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,71],type="l",pch=18,cex=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,75],type="l",pch=15,cex=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,79],type="l",pch=17,cex=1,col="orange",axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,83],type="l",pch=18,cex=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,87],type="l",pch=18,cex=1,col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,91],type="l",pch=18,cex=1,col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,95],type="l",pch=18,cex=1,col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " H3K4me3")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)

legend(90, 1.14, c("I","II","VI","VII","VIII"), lwd=1,
       col=c(rgb(0,0.7,1), rgb(0,0.2,0.5), rgb(0,0.6,0), rgb(1,0,1), rgb(0.55,0,0)),lty=c(1,1,6,6,6),
       pt.cex=1, y.intersp=0.45, x.intersp=0.1, seg.len=0.7, bty="n")
dev.off()


y[,1][which.max(y[,67])]
y[,1][which.max(y[,71])]
#y[,1][which.max(y[,75])]
#y[,1][which.max(y[,79])]
#y[,1][which.max(y[,83])]
y[,1][which.max(y[,87])]
y[,1][which.max(y[,91])]
y[,1][which.max(y[,95])]






# daily min
y <- read.csv("data/Histone_rep_SMAmin_120925-140916_log10sma_log10his_result.csv", header=T, sep=",")
y<-as.matrix(y)

# K27_STM
pdf("figs/R2_SMAmin-K27_STM_all.pdf", height=1.3, width=1.4)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,35],type="l",pch=16,cex=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,39],type="l",pch=18,cex=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,43],type="l",pch=15,cex=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,47],type="l",pch=17,cex=1,col="orange",axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,51],type="l",pch=18,cex=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,55],type="l",pch=18,cex=1,col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,59],type="l",pch=18,cex=1,col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,63],type="l",pch=18,cex=1,col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " H3K27me3")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)

legend(par()$usr[2]-35, par()$usr[4]+0.12, c("I","II","III","IV","V","VI","VII","VIII"), lwd=1,
       col=c(rgb(0,0.7,1), rgb(0,0.2,0.5), rgb(1,0,0.2), "orange", rgb(0.6,0.3,0.6), rgb(0,0.6,0), rgb(1,0,1), rgb(0.55,0,0)),lty=c(1,1,1,1,1,6,6,6),
       pt.cex=1, y.intersp=0.45, x.intersp=0.1, seg.len=0.7, bty="n")
dev.off()

y[,1][which.max(y[,35])]
y[,1][which.max(y[,39])]
y[,1][which.max(y[,43])]
y[,1][which.max(y[,47])]
y[,1][which.max(y[,51])]
y[,1][which.max(y[,55])]
y[,1][which.max(y[,59])]
y[,1][which.max(y[,63])]


# K4_ACT2
pdf("figs/R2_SMAmin-K4_ACT2_all.pdf", height=1.3, width=1.4)
par(mar=c(1.4,1.4,1,1.6))
par(mgp=c(0,0,0))
par(ps=6)
par(xpd=T)

plot(y[,1],y[,67],type="l",pch=16,cex=1,col=rgb(0,0.7,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,71],type="l",pch=18,cex=1,col=rgb(0,0.2,0.5),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,75],type="l",pch=15,cex=1,col=rgb(1,0,0.2),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,79],type="l",pch=17,cex=1,col="orange",axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
#par(new=T)
#plot(y[,1],y[,83],type="l",pch=18,cex=1,col=rgb(0.6,0.3,0.6),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1))
par(new=T)
plot(y[,1],y[,87],type="l",pch=18,cex=1,col=rgb(0,0.6,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,91],type="l",pch=18,cex=1,col=rgb(1,0,1),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)
par(new=T)
plot(y[,1],y[,95],type="l",pch=18,cex=1,col=rgb(0.55,0,0),axes=F,ann=F,xlim=c(0,200),ylim=c(0,1),lty=6)

axis(side=1, at=seq(0,200,40), lab=F, tcl=-0.1, las=1, pos=0)
mtext(c("0","40","80","120","160","200"), at=c(0,40,80,116,157,200), side=1, line=-0.5)
axis(side=2, at=seq(0,1,0.2), lab=F, tcl=-0.1, las=2, pos=0)
mtext(c("0.0","0.2","0.4","0.6","0.8","1.0"), at=seq(0,1,0.2), side=2, las=1, line=0)
mtext("Moving average", side=1, line=-0.07)
mtext("period (days)", side=1, line=0.35)
mtext(expression(paste(R^{2}, " value")), side=2, line=0.5)
mtext(expression(paste(italic("AhgFLC"), " H3K4me3")), side=3, line=-0.2)

arrows(0,1,200,1,length=0)
arrows(200,0,200,1,length=0)

legend(90, 1.14, c("I","II","VI","VII","VIII"), lwd=1,
       col=c(rgb(0,0.7,1), rgb(0,0.2,0.5), rgb(0,0.6,0), rgb(1,0,1), rgb(0.55,0,0)),lty=c(1,1,6,6,6),
       pt.cex=1, y.intersp=0.45, x.intersp=0.1, seg.len=0.7, bty="n")
dev.off()


y[,1][which.max(y[,67])]
y[,1][which.max(y[,71])]
#y[,1][which.max(y[,75])]
#y[,1][which.max(y[,79])]
#y[,1][which.max(y[,83])]
y[,1][which.max(y[,87])]
y[,1][which.max(y[,91])]
y[,1][which.max(y[,95])]

