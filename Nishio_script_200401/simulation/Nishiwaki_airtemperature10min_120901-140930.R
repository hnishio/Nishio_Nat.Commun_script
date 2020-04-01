
read.temperature <- function(csv) {
	data <- read.csv(csv)
	times <- as.POSIXct(paste(data[,'year'], data[,'month'], data[,'day'], data[,'time'], sep='/'), format='%Y/%m/%d/%H:%M')
	idx <- which(!is.na(times))
	idx2 <- idx[c(1, which(diff(times[idx])>0)+1)]
	idx3 <- idx2[!is.na(data[idx2, 'air.temperature'])]
	return(data.frame(time=times[idx3], temperature=data[idx3,'air.temperature']))
}
temperature <- read.temperature('data/180823_JMA10min_nishiwaki_120101-141231.csv')


temperature <- temperature[which(temperature$time=='2012-09-01 00:00:00'):which(temperature$time=='2014-10-01 00:00:00'),]

at.dates <- 
  c('2012-09-01','2012-10-01','2012-11-01','2012-12-01','2013-01-01','2013-02-01',
    '2013-03-01','2013-04-01','2013-05-01','2013-06-01','2013-07-01','2013-08-01',
    '2013-09-01','2013-10-01','2013-11-01','2013-12-01','2014-01-01','2014-02-01',
    '2014-03-01','2014-04-01','2014-05-01','2014-06-01','2014-07-01','2014-08-01',
    '2014-09-01','2014-10-01')




pdf("figs/Nishiwaki_airtemperature10min_120901-140930_190712.pdf",height=1.8,width=3.7)
	par(mar=c(1.7,2.1,1.2,1.2))
	par(mgp=c(0,0.1,0))
	par(ps=6)
	par(cex=1)
	par(xpd=T)
	
plot(NULL, axes=F,ann=F, 
    xlim=c(as.numeric(as.POSIXct("2012-09-01 00:00:00 JST")),
         as.numeric(as.POSIXct("2014-10-01 00:00:00 JST"))),
          xlab='', ylab='',ylim=c(-7,37))

rect(as.numeric(as.POSIXct("2012-09-01 00:00:00 JST")), 0, as.numeric(as.POSIXct("2014-10-01 00:00:00 JST")), 16, col=rgb(0.75,0.75,0.75,alpha=0.5), border="transparent")

par(new=T)
plot(temperature$time, temperature$temperature, 
    col="black", type='l', axes=F,ann=F, 
    xlim=c(as.numeric(as.POSIXct("2012-09-01 00:00:00 JST")),
         as.numeric(as.POSIXct("2014-10-01 00:00:00 JST"))), 
    lty=1, lwd=0.3, xlab='', ylab='',ylim=c(-7,37))
axis(side=2, at=seq(-5,35,5), las=2, tcl=0.2, pos=as.numeric(as.POSIXct("2012-09-01 JST")))
arrows(as.numeric(as.POSIXct("2012-09-01 JST")),-7,as.numeric(as.POSIXct("2012-09-01 JST")),37,length=0)

for(i in c(-5,0,5,10,15,20,25,30,35)){
	arrows(as.numeric(as.POSIXct("2012-09-01 JST")),i,as.numeric(as.POSIXct("2014-10-01 JST")),
	i,length=0,lwd=0.3)
}

for(i in c("2012-10-01 JST","2012-11-01 JST","2012-12-01 JST","2013-01-01 JST","2013-02-01 JST","2013-03-01 JST","2013-04-01 JST","2013-05-01 JST","2013-06-01 JST","2013-07-01 JST","2013-08-01 JST","2013-09-01 JST","2013-10-01 JST","2013-11-01 JST","2013-12-01 JST","2014-01-01 JST","2014-02-01 JST","2014-03-01 JST","2014-04-01 JST","2014-05-01 JST","2014-06-01 JST","2014-07-01 JST","2014-08-01 JST","2014-09-01 JST")){
	arrows(as.numeric(as.POSIXct(i)),-7,as.numeric(as.POSIXct(i)),
	37,length=0,lwd=0.3)
}

arrows(as.numeric(as.POSIXct("2012-09-01 00:00:00 JST")), 10.5, as.numeric(as.POSIXct("2014-10-01 00:00:00 JST")), 10.5,col="red",length=0,lwd=0.5)

mtext("Temperature (°C)",side=2,line=0.1)

place<-c(as.POSIXct(at.dates))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]),mean(place[13:14]),mean(place[14:15]),mean(place[15:16]),mean(place[16:17]),mean(place[17:18]),mean(place[18:19]),mean(place[19:20]),mean(place[20:21]),mean(place[21:22]),mean(place[22:23]),mean(place[23:24]),mean(place[24:25]),mean(place[25:26]))
labelx<-c("9","10","11","12","1","2","3","4","5","6","7","8","9","10","11","12","1","2","3","4","5","6","7","8","9")
axis(side=1, at=place, labels=F,tcl=-0.2, pos=-7)
mtext(side=1,at=lit,labelx,line=-0.5)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),37,as.numeric(as.POSIXct("2014-10-01 JST")),37,length=0)
arrows(as.numeric(as.POSIXct("2014-10-01 JST")),-7,as.numeric(as.POSIXct("2014-10-01 JST")),37,length=0)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),-12.5,
       as.numeric(as.POSIXct("2012-12-31 JST")),-12.5,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2013-01-01 JST")),-12.5,
       as.numeric(as.POSIXct("2013-12-31 JST")),-12.5,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2014-01-01 JST")),-12.5,
       as.numeric(as.POSIXct("2014-09-30 JST")),-12.5,length=0.03,code=3,lwd=1)
mtext("2012",side=1,at=as.numeric(as.POSIXct("2012-10-31 JST")),line=0.15)
mtext("2013",side=1,at=as.numeric(as.POSIXct("2013-06-30 JST")),line=0.15)
mtext("2014",side=1,at=as.numeric(as.POSIXct("2014-05-15 JST")),line=0.15)

dev.off()

