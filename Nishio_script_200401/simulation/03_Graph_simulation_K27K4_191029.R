
library(fields)
source('01_simulation_K27K4_191029.R')



# Read tempareature data
temperature <- read.temperature('data/180823_JMA10min_nishiwaki_120101-141231.csv')

# Read expression and histone modifications data
data <- read.csv('data/RNAACT2.K4ACT2.K27STM.190930.csv')

# Sampling time
dates <- as.POSIXct(
  paste(as.character(20000000+data[,'date']), "1200"), 
  format="%Y%m%d %H%M")

# Simulation start, end and step
time.start <- as.POSIXct('2012-09-25')
time.end <- dates[length(dates)]
time.step <- 60

# Trimming of temperature data
temperature <- 
  temperature[time.start<=temperature$time & 
                temperature$time<=time.end+as.difftime(1,units='days'), ]

# Function of simulation for given parameters
simulate.param <- function(
  model.fun, 
  parameter, 
  temperature, 
  time.start, time.step, time.end,
  init=NULL){
  
  # Temperature at start
  temperature.start <- temperature$temperature[temperature$time == time.start]
  
  # Setting of initial values
  initial.state <- c(uu=1, mu=0, mm=0, um=0, a1=1, a2=0)
  
  # Execute simulation
  sim.results <- simulate(
    initialState=initial.state,
    temperature=temperature,
    time.start=time.start,
    time.step=time.step,
    time.end=time.end,
    parameters=parameter,
    func=model.fun
  )
  
  return(sim.results)
}  


load('optim/optim_parameters_K27K4_model1_SANN_191029')
sim.model1 <- simulate.param(model1, as.list(parameter.optim), 
                            temperature, time.start, time.step, time.end)

load('optim/optim_parameters_K27K4_model2_SANN_191029')
sim.model2 <- simulate.param(model2, as.list(parameter.optim), 
                            temperature, time.start, time.step, time.end)

load('optim/optim_parameters_K27K4_model3_SANN_191029')
sim.model3 <- simulate.param(model3, as.list(parameter.optim), 
                            temperature, time.start, time.step, time.end)

load('optim/optim_parameters_K27K4_model4_SANN_191029')
sim.model4 <- simulate.param(model4, as.list(parameter.optim), 
                            temperature, time.start, time.step, time.end)




##### Visualization of the results #####

# Labels for x-axis
at.dates <- 
  c('2012-09-01','2012-10-01','2012-11-01','2012-12-01','2013-01-01','2013-02-01',
    '2013-03-01','2013-04-01','2013-05-01','2013-06-01','2013-07-01','2013-08-01',
    '2013-09-01','2013-10-01','2013-11-01','2013-12-01','2014-01-01','2014-02-01',
    '2014-03-01','2014-04-01','2014-05-01','2014-06-01','2014-07-01','2014-08-01',
    '2014-09-01','2014-10-01')

# Load result
m1.model1 <- max(data$K27.NR.STM)*(sim.model1$state[,'mu'] + sim.model1$state[,'mm'])
m2.model1 <- max(data$K27.distNR.STM)*(sim.model1$state[,'um'] + sim.model1$state[,'mm'])
a1.model1 <- max(data$K4.NR.ACT2)*(sim.model1$state[,'a1'])
a2.model1 <- max(data$K4.distNR.ACT2)*(sim.model1$state[,'a2'])

m1.model2 <- max(data$K27.NR.STM)*(sim.model2$state[,'mu'] + sim.model2$state[,'mm'])
m2.model2 <- max(data$K27.distNR.STM)*(sim.model2$state[,'um'] + sim.model2$state[,'mm'])
a1.model2 <- max(data$K4.NR.ACT2)*(sim.model2$state[,'a1'])
a2.model2 <- max(data$K4.distNR.ACT2)*(sim.model2$state[,'a2'])

m1.model3 <- max(data$K27.NR.STM)*(sim.model3$state[,'mu'] + sim.model3$state[,'mm'])
m2.model3 <- max(data$K27.distNR.STM)*(sim.model3$state[,'um'] + sim.model3$state[,'mm'])
a1.model3 <- max(data$K4.NR.ACT2)*(sim.model3$state[,'a1'])
a2.model3 <- max(data$K4.distNR.ACT2)*(sim.model3$state[,'a2'])

m1.model4 <- max(data$K27.NR.STM)*(sim.model4$state[,'mu'] + sim.model4$state[,'mm'])
m2.model4 <- max(data$K27.distNR.STM)*(sim.model4$state[,'um'] + sim.model4$state[,'mm'])
a1.model4 <- max(data$K4.NR.ACT2)*(sim.model4$state[,'a1'])
a2.model4 <- max(data$K4.distNR.ACT2)*(sim.model4$state[,'a2'])


# Read daily average of temperature
read.temperature <- function(csv) {
  data <- read.csv(csv)
  times <- as.POSIXct(paste(data[,'year'], data[,'month'], data[,'day'], data[,'time'], sep='/'), format='%Y/%m/%d/%H:%M')
  idx <- which(!is.na(times))
  idx <- idx[c(1, which(diff(times[idx])>0)+1)]
  idx <- idx[!is.na(data[idx, 'air.temperature'])]
  return(data.frame(time=times[idx], temperature=data[idx,'air.temperature']))
}
temperature <- read.temperature('data/120901-140930_dailymean_temp.csv')


place<-c(as.POSIXct(at.dates))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]),mean(place[13:14]),mean(place[14:15]),mean(place[15:16]),mean(place[16:17]),mean(place[17:18]),mean(place[18:19]),mean(place[19:20]),mean(place[20:21]),mean(place[21:22]),mean(place[22:23]),mean(place[23:24]),mean(place[24:25]),mean(place[25:26]))
labelx<-c("9","10","11","12","1","2","3","4","5","6","7","8","9","10","11","12","1","2","3","4","5","6","7","8","9")


dir.create("figs")

###### RNA expression

# log-log plot
plot(log10(data$K4.NR.ACT2), data$RNA.ACT2.log)

# linear regression
m <- lm(data$RNA.ACT2.log~log10(data$K4.NR.ACT2))

# plot results
x <- seq(0, 1.2, length.out=1000)
plot(data$K4.NR.ACT2, data$RNA.ACT2.log)
lines(x, m$coefficients[1] + m$coefficients[2]*log10(x))



pdf("figs/simulation_RNA_K27K4_model1_2_SANN_log10_scale_191029.pdf",height=1.2,width=3.7)
par(oma=c(0,0,0,0))
par(mar=c(1.7,2.1,1.2,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)
plot(temperature$time, temperature$temperature, 
     col="gray", type='l', axes=F,ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     lty=1, lwd=1, xlab='', ylab='',ylim=c(-5,35))
axis(side=4, at=seq(-5,35,10), las=2, tcl=0.2, pos=as.numeric(as.POSIXct("2014-10-01 JST")))
mtext("Temperature (°C)",side=4,line=-0.3)

par(xpd=T)
par(new=T)
plot(sim.model1$times, scale(m$coefficients[1] + m$coefficients[2]*log10(a1.model1)), 
	 type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(-4,2))

par(new=T)
plot(sim.model2$times, scale(m$coefficients[1] + m$coefficients[2]*log10(a1.model2)), 
	 type='l', lwd=1, las=2, tcl=0.3, col="red",
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(-4,2))

points(dates, scale(data$RNA.ACT2.log), lwd=1, col="black", pch=1,cex=0.7)

axis(side=1, at=place, labels=F,tcl=-0.2, pos=-4)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(-4,2,2),lab=c("-4","-2","0","2"),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
mtext("Standardised",side=2,line=0.45)
mtext(expression(paste(italic("AhgFLC")," mRNA",sep="")),side=2,line=-0.1)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),2,as.numeric(as.POSIXct("2014-09-30 JST")),2,length=0)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),-5.3,
       as.numeric(as.POSIXct("2012-12-31 JST")),-5.3,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2013-01-01 JST")),-5.3,
       as.numeric(as.POSIXct("2013-12-31 JST")),-5.3,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2014-01-01 JST")),-5.3,
       as.numeric(as.POSIXct("2014-09-30 JST")),-5.3,length=0.03,code=3,lwd=1)
mtext("2012",side=1,at=as.numeric(as.POSIXct("2012-10-31 JST")),line=0.15)
mtext("2013",side=1,at=as.numeric(as.POSIXct("2013-06-30 JST")),line=0.15)
mtext("2014",side=1,at=as.numeric(as.POSIXct("2014-05-15 JST")),line=0.15)

legend(as.numeric(as.POSIXct("2012-09-25 JST")),4.9,
       legend="Obs.",lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(as.numeric(as.POSIXct("2013-01-28 JST")),4.9,legend="Sim. (model 1, full)",
       lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(as.numeric(as.POSIXct("2013-10-20 JST")),4.9,legend=expression(paste("Sim. (model 2, without ",italic(a)[N],")",sep="")),
       col="red",lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
dev.off()



pdf("figs/simulation_RNA_K27K4_model1_3_SANN_log10_scale_191029.pdf",height=1.2,width=3.7)
par(oma=c(0,0,0,0))
par(mar=c(1.7,2.1,1.2,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)
plot(temperature$time, temperature$temperature, 
     col="gray", type='l', axes=F,ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     lty=1, lwd=1, xlab='', ylab='',ylim=c(-5,35))
axis(side=4, at=seq(-5,35,10), las=2, tcl=0.2, pos=as.numeric(as.POSIXct("2014-10-01 JST")))
mtext("Temperature (°C)",side=4,line=-0.3)

par(xpd=T)
par(new=T)
plot(sim.model1$times, scale(m$coefficients[1] + m$coefficients[2]*log10(a1.model1)), 
	 type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(-4,2))

par(new=T)
plot(sim.model3$times, scale(m$coefficients[1] + m$coefficients[2]*log10(a1.model3)), 
	 type='l', lwd=1, las=2, tcl=0.3, col="red",
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(-4,2))


points(dates, scale(data$RNA.ACT2.log), lwd=1, col="black", pch=1,cex=0.7)

axis(side=1, at=place, labels=F,tcl=-0.2, pos=-4)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(-4,2,2),lab=c("-4","-2","0","2"),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
mtext("Standardised",side=2,line=0.45)
mtext(expression(paste(italic("AhgFLC")," mRNA",sep="")),side=2,line=-0.1)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),2,as.numeric(as.POSIXct("2014-09-30 JST")),2,length=0)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),-5.3,
       as.numeric(as.POSIXct("2012-12-31 JST")),-5.3,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2013-01-01 JST")),-5.3,
       as.numeric(as.POSIXct("2013-12-31 JST")),-5.3,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2014-01-01 JST")),-5.3,
       as.numeric(as.POSIXct("2014-09-30 JST")),-5.3,length=0.03,code=3,lwd=1)
mtext("2012",side=1,at=as.numeric(as.POSIXct("2012-10-31 JST")),line=0.15)
mtext("2013",side=1,at=as.numeric(as.POSIXct("2013-06-30 JST")),line=0.15)
mtext("2014",side=1,at=as.numeric(as.POSIXct("2014-05-15 JST")),line=0.15)

legend(as.numeric(as.POSIXct("2012-08-11 JST")),4.9,
       legend="Obs.",lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(as.numeric(as.POSIXct("2012-11-20 JST")),4.9,legend="Sim. (model 1, full)",
       lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(as.numeric(as.POSIXct("2013-07-25 JST")),4.9,legend=expression(paste("Sim. [model 3, without (",italic(m)[N],italic(u)[D]," + ",italic(m)[N],italic(m)[D],")]",sep="")),
       col="red",lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
dev.off()



pdf("figs/simulation_RNA_K27K4_model1_4_SANN_log10_scale_191029.pdf",height=1.2,width=3.7)
par(oma=c(0,0,0,0))
par(mar=c(1.7,2.1,1.2,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)
plot(temperature$time, temperature$temperature, 
     col="gray", type='l', axes=F,ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     lty=1, lwd=1, xlab='', ylab='',ylim=c(-5,35))
axis(side=4, at=seq(-5,35,10), las=2, tcl=0.2, pos=as.numeric(as.POSIXct("2014-10-01 JST")))
mtext("Temperature (°C)",side=4,line=-0.3)

par(xpd=T)
par(new=T)
plot(sim.model1$times, scale(m$coefficients[1] + m$coefficients[2]*log10(a1.model1)), 
	 type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(-4,2))

par(new=T)
plot(sim.model4$times, scale(m$coefficients[1] + m$coefficients[2]*log10(a1.model4)), 
	 type='l', lwd=1, las=2, tcl=0.3, col="red",
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(-4,2))


points(dates, scale(data$RNA.ACT2.log), lwd=1, col="black", pch=1,cex=0.7)

axis(side=1, at=place, labels=F,tcl=-0.2, pos=-4)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(-4,2,2),lab=c("-4","-2","0","2"),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
mtext("Standardised",side=2,line=0.45)
mtext(expression(paste(italic("AhgFLC")," mRNA",sep="")),side=2,line=-0.1)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),2,as.numeric(as.POSIXct("2014-09-30 JST")),2,length=0)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),-5.3,
       as.numeric(as.POSIXct("2012-12-31 JST")),-5.3,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2013-01-01 JST")),-5.3,
       as.numeric(as.POSIXct("2013-12-31 JST")),-5.3,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2014-01-01 JST")),-5.3,
       as.numeric(as.POSIXct("2014-09-30 JST")),-5.3,length=0.03,code=3,lwd=1)
mtext("2012",side=1,at=as.numeric(as.POSIXct("2012-10-31 JST")),line=0.15)
mtext("2013",side=1,at=as.numeric(as.POSIXct("2013-06-30 JST")),line=0.15)
mtext("2014",side=1,at=as.numeric(as.POSIXct("2014-05-15 JST")),line=0.15)

legend(as.numeric(as.POSIXct("2012-08-11 JST")),4.9,
       legend="Obs.",lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(as.numeric(as.POSIXct("2012-11-20 JST")),4.9,legend="Sim. (model 1, full)",
       lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(as.numeric(as.POSIXct("2013-07-25 JST")),4.9,legend=expression(paste("Sim. [model 4, without (",italic(m)[N],italic(m)[D]," + ",italic(u)[N],italic(m)[D],")]",sep="")),
       col="red",lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
dev.off()




# Histone state

#model1
uu.model1 <- sim.model1$state[,'uu']
mu.model1 <- sim.model1$state[,'mu']
mm.model1 <- sim.model1$state[,'mm']
um.model1 <- sim.model1$state[,'um']

pdf("figs/simulation_histone_state_K27K4_model1_SANN_191029.pdf",height=1.2,width=3.7)
par(oma=c(0,0,0,0))
par(mar=c(1.7,2.1,1.2,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)

plot(temperature$time, temperature$temperature, 
     col="gray", type='l', axes=F,ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     lty=1, lwd=1, xlab='', ylab='',ylim=c(-5,35))
axis(side=4, at=seq(-5,35,10), las=2, tcl=0.2, pos=as.numeric(as.POSIXct("2014-10-01 JST")))
mtext("Temperature (°C)",side=4,line=-0.3)

par(xpd=T)
par(new=T)
plot(sim.model1$times, uu.model1, type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),  
     ylim=c(0,1.2))
par(new=T)
plot(sim.model1$times, mu.model1, type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),  
     ylim=c(0,1.2), col="red")
par(new=T)
plot(sim.model1$times, mm.model1, type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),  
     ylim=c(0,1.2), col="blue")
par(new=T)
plot(sim.model1$times, um.model1, type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),  
     ylim=c(0,1.2), col="magenta")

axis(side=1, at=place, labels=F,tcl=-0.2, pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1.2,0.4),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
#axis(side=2,at=seq(-0.1,1,0.1),lab=c("","0","","0.2","","0.4","","0.6","","0.8","","1"),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
mtext("Proportion (model 1)",side=2,line=0.15)
arrows(as.numeric(as.POSIXct("2012-09-01 JST")),1.2,as.numeric(as.POSIXct("2014-09-30 JST")),1.2,length=0)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),-0.26,
       as.numeric(as.POSIXct("2012-12-31 JST")),-0.26,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2013-01-01 JST")),-0.26,
       as.numeric(as.POSIXct("2013-12-31 JST")),-0.26,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2014-01-01 JST")),-0.26,
       as.numeric(as.POSIXct("2014-09-30 JST")),-0.26,length=0.03,code=3,lwd=1)
mtext("2012",side=1,at=as.numeric(as.POSIXct("2012-10-31 JST")),line=0.15)
mtext("2013",side=1,at=as.numeric(as.POSIXct("2013-06-30 JST")),line=0.15)
mtext("2014",side=1,at=as.numeric(as.POSIXct("2014-05-15 JST")),line=0.15)

legend(as.numeric(as.POSIXct("2013-01-17 JST")),1.75,
	legend=expression(paste(U[N],U[D],sep="")),lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",
	x.intersp=0.1,seg.len=0.7)
legend(as.numeric(as.POSIXct("2013-05-15 JST")),1.75,
	legend=expression(paste(M[N],U[D],sep="")),col="red",lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",
	x.intersp=0.1,seg.len=0.7)
legend(as.numeric(as.POSIXct("2013-09-13 JST")),1.75,
    legend=expression(paste(M[N],M[D],sep="")),col="blue",lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",
    x.intersp=0.1,seg.len=0.7)
legend(as.numeric(as.POSIXct("2014-01-13 JST")),1.75,
	legend=expression(paste(U[N],M[D],sep="")),col="magenta",lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",
	x.intersp=0.1,seg.len=0.7)

dev.off()




##### model1

# H3K27me3 at NR and distNR
pdf("figs/simulation_K27_K27K4_model1_SANN_191029-2.pdf",height=1.2,width=3.7)
par(oma=c(0,0,0,0))
par(mar=c(1.7,2.1,1.2,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)

plot(temperature$time, temperature$temperature, 
     col="gray", type='l', axes=F,ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     lty=1, lwd=1, xlab='', ylab='',ylim=c(-5,35))
axis(side=4, at=seq(-5,35,10), las=2, tcl=0.2, pos=as.numeric(as.POSIXct("2014-10-01 JST")))
mtext("Temperature (°C)", side = 4, line = -0.3)

par(xpd=T)
par(new=T)
plot(sim.model1$times, m1.model1, type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     ylim=c(0,1), col=rgb(0,0.7,1))
par(new=T)
plot(sim.model1$times, m2.model1, type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(0,1), col=rgb(0.6,0.3,0.6))

points(dates, data$K27.NR.STM, lwd=1, col=rgb(0,0.7,1), pch=1, cex=0.7)
points(dates, data$K27.distNR.STM, lwd=1, col=rgb(0.6,0.3,0.6), pch=1, cex=0.7)

axis(side=1, at=place, labels=F,tcl=-0.2, pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1,0.2),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
#axis(side=2,at=seq(-0.1,1,0.1),lab=c("","0","","0.2","","0.4","","0.6","","0.8","","1"),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.5)
mtext("H3K27me3",side=2,line=0.1)
arrows(as.numeric(as.POSIXct("2012-09-01 JST")),1,as.numeric(as.POSIXct("2014-09-30 JST")),1,length=0)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),-0.22,
       as.numeric(as.POSIXct("2012-12-31 JST")),-0.22,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2013-01-01 JST")),-0.22,
       as.numeric(as.POSIXct("2013-12-31 JST")),-0.22,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2014-01-01 JST")),-0.22,
       as.numeric(as.POSIXct("2014-09-30 JST")),-0.22,length=0.03,code=3,lwd=1)
mtext("2012",side=1,at=as.numeric(as.POSIXct("2012-10-31 JST")),line=0.15)
mtext("2013",side=1,at=as.numeric(as.POSIXct("2013-06-30 JST")),line=0.15)
mtext("2014",side=1,at=as.numeric(as.POSIXct("2014-05-15 JST")),line=0.15)

legend(as.numeric(as.POSIXct("2012-05-20 JST")),1.53,legend="Obs. (NR)",
       col=rgb(0,0.7,1),lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(as.numeric(as.POSIXct("2012-10-27 JST")),1.53,legend="Sim. (model 1, NR)",
       col=rgb(0,0.7,1),lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(as.numeric(as.POSIXct("2013-07-01 JST")),1.53,legend="Obs. (distal NR)",
       col=rgb(0.6,0.3,0.6),lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(as.numeric(as.POSIXct("2014-02-04 JST")),1.53,legend="Sim. (model 1, distal NR)",
       col=rgb(0.6,0.3,0.6),lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
dev.off()




# H3K4me3 at NR and distNR
pdf("figs/simulation_K4_K27K4_model1_SANN_191029-2.pdf",height=1.2,width=3.7)
par(oma=c(0,0,0,0))
par(mar=c(1.7,2.1,1.2,1.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(cex=1)

plot(temperature$time, temperature$temperature, 
     col="gray", type='l', axes=F,ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     lty=1, lwd=1, xlab='', ylab='',ylim=c(-5,35))
axis(side=4, at=seq(-5,35,10), las=2, tcl=0.2, pos=as.numeric(as.POSIXct("2014-10-01 JST")))
mtext("Temperature (°C)", side = 4, line = -0.3)

par(xpd=T)
par(new=T)
plot(sim.model1$times, a1.model1, type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     ylim=c(0,1.2), col=rgb(0,0.7,1))

par(new=T)
plot(sim.model1$times, a2.model1, type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     ylim=c(0,1.2), col=rgb(0.6,0.3,0.6))

points(dates, data$K4.NR.ACT2, lwd=1, col=rgb(0,0.7,1), pch=1, cex=0.7)
points(dates, data$K4.distNR.ACT2, lwd=1, col=rgb(0.6,0.3,0.6), pch=1, cex=0.7)

axis(side=1, at=place, labels=F,tcl=-0.2, pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,1.2,0.4),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
#axis(side=2,at=seq(-0.1,1,0.1),lab=c("","0","","0.2","","0.4","","0.6","","0.8","","1"),tcl=0.2,las=2,pos=as.numeric(as.POSIXct("2012-09-01 JST")))
mtext(expression(paste(italic("AhgFLC"))),side=2,line=0.5)
mtext("H3K4me3",side=2,line=0.1)
arrows(as.numeric(as.POSIXct("2012-09-01 JST")),1.2,as.numeric(as.POSIXct("2014-09-30 JST")),1.2,length=0)

arrows(as.numeric(as.POSIXct("2012-09-01 JST")),-0.264,
       as.numeric(as.POSIXct("2012-12-31 JST")),-0.264,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2013-01-01 JST")),-0.264,
       as.numeric(as.POSIXct("2013-12-31 JST")),-0.264,length=0.03,code=3,lwd=1)
arrows(as.numeric(as.POSIXct("2014-01-01 JST")),-0.264,
       as.numeric(as.POSIXct("2014-09-30 JST")),-0.264,length=0.03,code=3,lwd=1)
mtext("2012",side=1,at=as.numeric(as.POSIXct("2012-10-31 JST")),line=0.15)
mtext("2013",side=1,at=as.numeric(as.POSIXct("2013-06-30 JST")),line=0.15)
mtext("2014",side=1,at=as.numeric(as.POSIXct("2014-05-15 JST")),line=0.15)

legend(as.numeric(as.POSIXct("2012-05-20 JST")),1.836,legend="Obs. (NR)",
       col=rgb(0,0.7,1),lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(as.numeric(as.POSIXct("2012-10-27 JST")),1.836,legend="Sim. (model 1, NR)",
       col=rgb(0,0.7,1),lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(as.numeric(as.POSIXct("2013-07-01 JST")),1.836,legend="Obs. (distal NR)",
       col=rgb(0.6,0.3,0.6),lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(as.numeric(as.POSIXct("2014-02-04 JST")),1.836,legend="Sim. (model 1, distal NR)",
       col=rgb(0.6,0.3,0.6),lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
       
dev.off()
