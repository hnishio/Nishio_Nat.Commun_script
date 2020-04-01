
##### Calculation of Simple Moving Average #####
library(TTR)


y<-read.csv("data/2012-2014_temp.csv",header=T,sep=",")
head(y)

# SMA at sampling date
temp.ave<-matrix(ncol=200,nrow=50)
samplepoint <- c(269,297,311,326,339,353,374,388,402,416,430,444,458,472,
				 486,499,514,528,549,563,577,591,605,619,633,647,661,675,
				 689,703,717,738,752,766,780,794,810,822,836,849,864,878,
				 892,905,920,934,948,962,976,990)
samplepoint2 <- samplepoint-1

length(samplepoint2)

# daily mean
for(i in 1:200){
 	mvtemp<-SMA(y$daily_mean_temp,i)
 	for(j in 1:50){
 		temp.ave[j,i] <- mvtemp[samplepoint2[j]]
 		}
 	}
colnames(temp.ave) <- paste0("SMA",1:200)
write.csv(temp.ave,"data/SMAmean_120925-140916.csv",quote=F,row.names=F)

# daily max
temp.ave<-matrix(ncol=200,nrow=50)
for(i in 1:200){
  mvtemp<-SMA(y$daily_max_temp,i)
  for(j in 1:50){
    temp.ave[j,i] <- mvtemp[samplepoint2[j]]
  }
}
colnames(temp.ave) <- paste0("SMA",1:200)
write.csv(temp.ave,"data/SMAmax_120925-140916.csv",quote=F,row.names=F)

# daily min
temp.ave<-matrix(ncol=200,nrow=50)
for(i in 1:200){
  mvtemp<-SMA(y$daily_min_temp,i)
  for(j in 1:50){
    temp.ave[j,i] <- mvtemp[samplepoint2[j]]
  }
}
colnames(temp.ave) <- paste0("SMA",1:200)
write.csv(temp.ave,"data/SMAmin_120925-140916.csv",quote=F,row.names=F)




# SMA for visualisation (120901-140930)
temp.ave<-matrix(ncol=200,nrow=760)
samplepoint <- c(245:1004)
samplepoint2<-samplepoint-1

length(samplepoint2)

for(i in 1:200){
 	mvtemp<-SMA(y$temp,i)
 	for(j in 1:760){
 		temp.ave[j,i] <- mvtemp[samplepoint2[j]]
 		}
 	}
write.table(temp.ave,"data/SMA_120901-140930.csv",quote=F,col.names=F,row.names=F,sep=",")



####### SMA plot ################

x<-read.csv("data/SMA_120901-140930.csv",header=F,sep=",")
x2<-as.matrix(x)

#1d,1w,6w,12w,24w
j <- c(1,7,42,84,168)

tempcol <- c("orange","blue","red",rgb(0.4,0.8,1),"gray")

pdf("figs/SMA_temperature.pdf", height=1.2, width=3.7)
par(oma=c(0,0,0,0))
par(mar=c(1.7,2.1,1.2,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(xpd=T)

plot(1:760,x2[,j[1]],type="l",pch=0,col=tempcol[1],ann=F,axes=F,xlim=c(0,760),ylim=c(-5,35),cex=0.7)
par(new=T)
plot(1:760,x2[,j[2]],type="l",pch=15,col=tempcol[2],ann=F,axes=F,xlim=c(0,760),ylim=c(-5,35),cex=0.7)
par(new=T)
plot(1:760,x2[,j[3]],type="l",pch=2,col=tempcol[3],ann=F,axes=F,xlim=c(0,760),ylim=c(-5,35),cex=0.7)
par(new=T)
plot(1:760,x2[,j[4]],type="l",pch=17,col=tempcol[4],ann=F,axes=F,xlim=c(0,760),ylim=c(-5,35),cex=0.7)
par(new=T)
plot(1:760,x2[,j[5]],type="l",pch=4,col=tempcol[5],ann=F,axes=F,xlim=c(0,760),ylim=c(-5,35),cex=0.7)

day<-c(0,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),sum(day[1:11]),sum(day[1:12]),sum(day[1:13]),sum(day[1:14]),sum(day[1:15]),sum(day[1:16]),sum(day[1:17]),sum(day[1:18]),sum(day[1:19]),sum(day[1:20]),sum(day[1:21]),sum(day[1:22]),sum(day[1:23]),sum(day[1:24]),sum(day[1:25]),sum(day[1:26]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]),mean(place[13:14]),mean(place[14:15]),mean(place[15:16]),mean(place[16:17]),mean(place[17:18]),mean(place[18:19]),mean(place[19:20]),mean(place[20:21]),mean(place[21:22]),mean(place[22:23]),mean(place[23:24]),mean(place[24:25]),mean(place[25:26]))
labelx<-c("9","10","11","12","1","2","3","4","5","6","7","8","9","10","11","12","1","2","3","4","5","6","7","8","9")
axis(side=1,at=place,label=F,tcl=-0.2,pos=-5)
mtext(side=1,at=lit,labelx,line=-0.4)
axis(side=2,at=seq(-5,35,10),tcl=0.2,las=2,pos=0)
mtext("SMA of",side=2,line=0.5)
mtext("temperature (°C)",side=2,line=0.08)
arrows(0,35,760,35,length=0)
arrows(760,-5,760,35,length=0)
legend(-20, 55, "1 day", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", col = tempcol[1], pt.cex=0.7)
legend(103, 55, "1 week", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", col = tempcol[2], pt.cex=0.7)
legend(245, 55, "6 weeks", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", col = tempcol[3], pt.cex=0.7)
legend(400, 55, "12 weeks", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", col = tempcol[4], pt.cex=0.7)
legend(565, 55, "24 weeks", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", col = tempcol[5], pt.cex=0.7)

arrows(0,-15,122,-15,length=0.05,code=3,lwd=0.5)
arrows(122,-15,487,-15,length=0.05,code=3,lwd=0.5)
arrows(487,-15,760,-15,length=0.05,code=3,lwd=0.5)
mtext("2012",side=1,at=61,line=0.3)
mtext("2013",side=1,at=303,line=0.3)
mtext("2014",side=1,at=622,line=0.3)

dev.off()
