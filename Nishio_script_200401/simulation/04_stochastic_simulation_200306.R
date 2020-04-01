
set.seed(123)

##### Definition of stochastic models #####
# Transition to next state stochastically
# K27 UU:0, MU:1, MM:2, UM:3
# K4/K4D U:0, A:1
model1 <- function(K27, K4, K4D, temperature, parameters, N){
  trans.prob.K27.base <- c(
    parameters$zeta / (1 + exp(parameters$alpha * (temperature - 5))), 
    parameters$eta / (1 + exp(-parameters$beta * (temperature - 10))), 
    parameters$kappa,
    parameters$lambda
  )
  idx.mm <- which(K27 == 2)
  trans.prob.K27 <- trans.prob.K27.base[K27 + 1]
  trans.prob.K27[idx.mm] <- trans.prob.K27[idx.mm] * K4[idx.mm]
  
  is.a1 <- (K4==1)
  is.a2 <- (K4D==1)
  trans.prob.K4 <- rep(0, N)
  trans.prob.K4D <- rep(0, N)
  trans.prob.K4[is.a1] <- parameters$phi * (K27[is.a1]==1 | K27[is.a1]==2)
  trans.prob.K4[!is.a1] <- parameters$rho / (1 + exp(-parameters$epsilon * (temperature - 15))) * (!is.a2[!is.a1])
  
  trans.prob.K4D[is.a2] <- parameters$psi * (K27[is.a2]==2 | K27[is.a2]==3)
  trans.prob.K4D[!is.a2] <- parameters$iota / (1 + exp(parameters$gamma * (temperature - 5))) * (!is.a1[!is.a2])
  
  list(
    K27 = (K27 + (runif(N)<trans.prob.K27)) %% 4, 
    K4  = (K4 + (runif(N)<trans.prob.K4)) %% 2, 
    K4D = (K4D + (runif(N)<trans.prob.K4D)) %%2
  )
}

model2 <- function(K27, K4, K4D, temperature, parameters, N){
  trans.prob.K27.base <- c(
    parameters$zeta / (1 + exp(parameters$alpha * (temperature - 5))), 
    parameters$eta / (1 + exp(-parameters$beta * (temperature - 10))), 
    parameters$kappa,
    parameters$lambda
  )
  trans.prob.K27 <- trans.prob.K27.base[K27 + 1]

  is.a1 <- (K4==1)
  is.a2 <- (K4D==1)
  trans.prob.K4 <- rep(0, N)
  trans.prob.K4D <- rep(0, N)
  trans.prob.K4[is.a1] <- parameters$phi * (K27[is.a1]==1 | K27[is.a1]==2)
  trans.prob.K4[!is.a1] <- parameters$rho / (1 + exp(-parameters$epsilon * (temperature - 15))) * (!is.a2[!is.a1])
  
  trans.prob.K4D[is.a2] <- parameters$psi * (K27[is.a2]==2 | K27[is.a2]==3)
  trans.prob.K4D[!is.a2] <- parameters$iota / (1 + exp(parameters$gamma * (temperature - 5))) * (!is.a1[!is.a2])
  
  list(
    K27 = (K27 + (runif(N)<trans.prob.K27)) %% 4, 
    K4  = (K4 + (runif(N)<trans.prob.K4)) %% 2, 
    K4D = (K4D + (runif(N)<trans.prob.K4D)) %%2
  )
}

model3 <- function(K27, K4, K4D, temperature, parameters, N){
  trans.prob.K27.base <- c(
    parameters$zeta / (1 + exp(parameters$alpha * (temperature - 5))), 
    parameters$eta / (1 + exp(-parameters$beta * (temperature - 10))), 
    parameters$kappa,
    parameters$lambda
  )
  idx.mm <- which(K27 == 2)
  trans.prob.K27 <- trans.prob.K27.base[K27 + 1]
  trans.prob.K27[idx.mm] <- trans.prob.K27[idx.mm] * K4[idx.mm]
  
  is.a1 <- (K4==1)
  is.a2 <- (K4D==1)
  trans.prob.K4 <- rep(0, N)
  trans.prob.K4D <- rep(0, N)
  trans.prob.K4[is.a1] <- parameters$phi
  trans.prob.K4[!is.a1] <- parameters$rho / (1 + exp(-parameters$epsilon * (temperature - 15))) * (!is.a2[!is.a1])
  
  trans.prob.K4D[is.a2] <- parameters$psi * (K27[is.a2]==2 | K27[is.a2]==3)
  trans.prob.K4D[!is.a2] <- parameters$iota / (1 + exp(parameters$gamma * (temperature - 5))) * (!is.a1[!is.a2])
  
  list(
    K27 = (K27 + (runif(N)<trans.prob.K27)) %% 4, 
    K4  = (K4 + (runif(N)<trans.prob.K4)) %% 2, 
    K4D = (K4D + (runif(N)<trans.prob.K4D)) %%2
  )
}
model4 <- function(K27, K4, K4D, temperature, parameters, N){
  trans.prob.K27.base <- c(
    parameters$zeta / (1 + exp(parameters$alpha * (temperature - 5))), 
    parameters$eta / (1 + exp(-parameters$beta * (temperature - 10))), 
    parameters$kappa,
    parameters$lambda
  )
  idx.mm <- which(K27 == 2)
  trans.prob.K27 <- trans.prob.K27.base[K27 + 1]
  trans.prob.K27[idx.mm] <- trans.prob.K27[idx.mm] * K4[idx.mm]
  
  is.a1 <- (K4==1)
  is.a2 <- (K4D==1)
  trans.prob.K4 <- rep(0, N)
  trans.prob.K4D <- rep(0, N)
  trans.prob.K4[is.a1] <- parameters$phi * (K27[is.a1]==1 | K27[is.a1]==2)
  trans.prob.K4[!is.a1] <- parameters$rho / (1 + exp(-parameters$epsilon * (temperature - 15))) * (!is.a2[!is.a1])
  
  trans.prob.K4D[is.a2] <- parameters$psi
  trans.prob.K4D[!is.a2] <- parameters$iota / (1 + exp(parameters$gamma * (temperature - 5))) * (!is.a1[!is.a2])
  
  list(
    K27 = (K27 + (runif(N)<trans.prob.K27)) %% 4, 
    K4  = (K4 + (runif(N)<trans.prob.K4)) %% 2, 
    K4D = (K4D + (runif(N)<trans.prob.K4D)) %%2
  )
}

# Read the parameters optimised in ODE models
parameters <- list()
for(i in 1:4){
  load(sprintf("optim/optim_parameters_K27K4_model%d_SANN_191029", i))
  # Adjust time interval to set one step to one hour
  parameters[[i]] <-as.list(parameter.optim)
  parameters[[i]]$zeta <- parameters[[i]]$zeta / 1440
  parameters[[i]]$eta <- parameters[[i]]$eta / 1440
  parameters[[i]]$iota <- parameters[[i]]$iota / 1440
  parameters[[i]]$rho <- parameters[[i]]$rho / 1440
  parameters[[i]]$lambda <- parameters[[i]]$lambda / 1440
  parameters[[i]]$kappa <- parameters[[i]]$kappa / 1440
  parameters[[i]]$phi <- parameters[[i]]$phi / 1440
  parameters[[i]]$psi <- parameters[[i]]$psi / 1440
}

# Load the function
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
times <- seq(time.start, time.end, by=time.step*60)

# Index for sampling time
data.idx <- as.numeric((dates - time.start) / time.step, units='mins')

# Trimming of temperature data
temperature <- temperature[time.start<=temperature$time & 
                             temperature$time<=time.end+as.difftime(1,units='days'), ]
temp.fun <- approxfun(as.numeric(temperature$time - time.start), temperature$temperature)
temp <- temp.fun(times - times[1])

# Simulate
models <- list(model1, model2, model3, model4)
N <- 10000
L <- length(times)
s.state <- list()
s.raw <- list()
for(i in 1:4){
  s.raw[[i]] <- list()
  s.raw[[i]]$K27 <- rep(0, N)
  s.raw[[i]]$K4 <- rep(1, N)
  s.raw[[i]]$K4D <- rep(0, N)
  s.state[[i]] <- matrix(c(1, 0, 0, 0, 1, 0), L, 6, byrow = T)
  colnames(s.state[[i]]) <- c("uu", "mu", "mm", "um", "a1", "a2")
}
for(i in 2:L){
  if(i %% 100 == 0) cat(sprintf("%d/%d\n", i, L))
  for(j in 1:4){
    s <- models[[j]](s.raw[[j]]$K27, s.raw[[j]]$K4, s.raw[[j]]$K4D, temp[i-1], parameters[[j]], N)
    s.raw[[j]]$K27 <- s$K27
    s.raw[[j]]$K4 <- s$K4
    s.raw[[j]]$K4D <- s$K4D
    s.state[[j]][i,] <- c(sapply(0:3, function(i)mean(s.raw[[j]]$K27==i)), mean(s.raw[[j]]$K4), mean(s.raw[[j]]$K4D))
  }
}
# Save the results
save(s.state, file="stochastic_results_raw")


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

# Load the results of ODE models
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

m <- lm(data$RNA.ACT2.log~log10(data$K4.NR.ACT2))

load("stochastic_results_raw")



##### Visualization of the results #####
sim.models <- rbind(
data.frame(time=times, sim.model1$state[,-1], simulation="ODE", model="model1"),
data.frame(time=times, sim.model2$state[,-1],simulation="ODE", model="model2"),
data.frame(time=times, sim.model3$state[,-1],simulation="ODE", model="model3"),
data.frame(time=times, sim.model4$state[,-1],simulation="ODE", model="model4"),
data.frame(time=times, s.state[[1]],simulation="Stochastic", model="model1"),
data.frame(time=times, s.state[[2]],simulation="Stochastic", model="model2"),
data.frame(time=times, s.state[[3]],simulation="Stochastic", model="model3"),
data.frame(time=times, s.state[[4]],simulation="Stochastic", model="model4")
)

sims <- data.frame(sim.models, 
           mRNA=m$coefficients[1] + m$coefficients[2]*log10(sim.models$a1),
           m1=max(data$K27.NR.STM)*(sim.models$mu+sim.models$mm), 
           m2=max(data$K27.distNR.STM)*(sim.models$um+sim.models$mm),
           a1.rev=max(data$K4.NR.ACT2)*sim.models$a1, 
           a2.rev=max(data$K4.distNR.ACT2)*sim.models$a2)
sims2 <- data.frame(sims,
                   scaled.mRNA=unlist(tapply(sims$mRNA, 
                                             paste(sims$simulation, sims$model, sep="-"), 
                                             scale)
                                      )
                   )

data2 <- data.frame(data, time=dates, scaled.mRNA=scale(data$RNA.ACT2.log))


# Labels for x-axis
at.dates <- 
  c('2012-09-01','2012-10-01','2012-11-01','2012-12-01','2013-01-01','2013-02-01',
    '2013-03-01','2013-04-01','2013-05-01','2013-06-01','2013-07-01','2013-08-01',
    '2013-09-01','2013-10-01','2013-11-01','2013-12-01','2014-01-01','2014-02-01',
    '2014-03-01','2014-04-01','2014-05-01','2014-06-01','2014-07-01','2014-08-01',
    '2014-09-01','2014-10-01')

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


# mRNA
pdf("figs/stochastic_RNA_K27K4_model1_2_SANN_log10_scale_200304.pdf",height=1.2,width=3.7)
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
plot(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
     sims2$scaled.mRNA[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
     type='l', lwd=1, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(-4,2))

par(new=T)
plot(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model2"],
     sims2$scaled.mRNA[sims2$simulation=="Stochastic"&sims2$model=="model2"], 
     type='l', lwd=1, las=2, tcl=0.3, col="red",
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(-4,2))

points(data2$time, data2$scaled.mRNA, lwd=1, col="black", pch=1,cex=0.7)

place<-c(as.POSIXct(at.dates))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]),mean(place[13:14]),mean(place[14:15]),mean(place[15:16]),mean(place[16:17]),mean(place[17:18]),mean(place[18:19]),mean(place[19:20]),mean(place[20:21]),mean(place[21:22]),mean(place[22:23]),mean(place[23:24]),mean(place[24:25]),mean(place[25:26]))
labelx<-c("9","10","11","12","1","2","3","4","5","6","7","8","9","10","11","12","1","2","3","4","5","6","7","8","9")
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

#legend(as.numeric(as.POSIXct("2012-09-25 JST")),4.9,
#       legend="Obs.",lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
#legend(as.numeric(as.POSIXct("2013-01-28 JST")),4.9,legend="Sim. (model 1, full)",
#       lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
#legend(as.numeric(as.POSIXct("2013-10-20 JST")),4.9,legend=expression(paste("Sim. (model 2, without ",italic(a)[N],")",sep="")),
#       col="red",lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
dev.off()


# Histone state
#model1
pdf("figs/stochastic_histone_state_K27K4_model1_SANN_200304.pdf",height=1.2,width=3.7)
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
plot(sims2$time[sims2$simulation=="ODE"&sims2$model=="model1"], 
     sims2$uu[sims2$simulation=="ODE"&sims2$model=="model1"], 
     type='l', lwd=0.5, las=2,tcl=0.3, xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),  
     ylim=c(0,1.2))
lines(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
      sims2$uu[sims2$simulation=="Stochastic"&sims2$model=="model1"],
      type="l", lwd=1.5)
par(new=T)
plot(sims2$time[sims2$simulation=="ODE"&sims2$model=="model1"],
     sims2$mu[sims2$simulation=="ODE"&sims2$model=="model1"], 
     type='l', lwd=0.5, las=2,tcl=0.3, xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),  
     ylim=c(0,1.2), col="red")
lines(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
      sims2$mu[sims2$simulation=="Stochastic"&sims2$model=="model1"],
      type="l", col="red", lwd=1.5)
par(new=T)
plot(sims2$time[sims2$simulation=="ODE"&sims2$model=="model1"],
     sims2$mm[sims2$simulation=="ODE"&sims2$model=="model1"], 
     type='l', lwd=0.5, las=2,tcl=0.3, xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),  
     ylim=c(0,1.2), col="blue")
lines(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
      sims2$mm[sims2$simulation=="Stochastic"&sims2$model=="model1"],
      type="l", col="blue", lwd=1.5)
par(new=T)
plot(sims2$time[sims2$simulation=="ODE"&sims2$model=="model1"],
     sims2$um[sims2$simulation=="ODE"&sims2$model=="model1"], 
     type='l', lwd=0.5, las=2,tcl=0.3, xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),  
     ylim=c(0,1.2), col="magenta")
lines(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
      sims2$um[sims2$simulation=="Stochastic"&sims2$model=="model1"],
      type="l", col="magenta", lwd=1.5)

place<-c(as.POSIXct(at.dates))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]),mean(place[13:14]),mean(place[14:15]),mean(place[15:16]),mean(place[16:17]),mean(place[17:18]),mean(place[18:19]),mean(place[19:20]),mean(place[20:21]),mean(place[21:22]),mean(place[22:23]),mean(place[23:24]),mean(place[24:25]),mean(place[25:26]))
labelx<-c("9","10","11","12","1","2","3","4","5","6","7","8","9","10","11","12","1","2","3","4","5","6","7","8","9")
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
dev.off()


# H3K27me3 at NR and distNR
pdf("figs/stochastic_K27_K27K4_model1_SANN_200304.pdf",height=1.2,width=3.7)
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
plot(sims2$time[sims2$simulation=="ODE"&sims2$model=="model1"], 
     sims2$m1[sims2$simulation=="ODE"&sims2$model=="model1"],
     type='l', lwd=0.5, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     ylim=c(0,1), col=rgb(0,0.7,1))
lines(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
      sims2$m1[sims2$simulation=="Stochastic"&sims2$model=="model1"],
      type="l", lwd=1.5, col=rgb(0,0.7,1))
par(new=T)
plot(sims2$time[sims2$simulation=="ODE"&sims2$model=="model1"], 
     sims2$m2[sims2$simulation=="ODE"&sims2$model=="model1"],
     type='l', lwd=0.5, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))),
     ylim=c(0,1), col=rgb(0.6,0.3,0.6))
lines(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
      sims2$m2[sims2$simulation=="Stochastic"&sims2$model=="model1"],
      type="l", lwd=1.5, col=rgb(0.6,0.3,0.6))

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
dev.off()


# H3K4me3 at NR and distNR
pdf("figs/stochastic_K4_K27K4_model1_SANN_200304.pdf",height=1.2,width=3.7)
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
plot(sims2$time[sims2$simulation=="ODE"&sims2$model=="model1"], 
     sims2$a1.rev[sims2$simulation=="ODE"&sims2$model=="model1"],
     type='l', lwd=0.5, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     ylim=c(0,1.2), col=rgb(0,0.7,1))
lines(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
      sims2$a1.rev[sims2$simulation=="Stochastic"&sims2$model=="model1"],
      type="l", lwd=1.5, col=rgb(0,0.7,1))

par(new=T)
plot(sims2$time[sims2$simulation=="ODE"&sims2$model=="model1"], 
     sims2$a2.rev[sims2$simulation=="ODE"&sims2$model=="model1"],
     type='l', lwd=0.5, las=2,tcl=0.3,
     xaxt='n', axes=F, ann=F, 
     xlim=c(as.numeric(as.POSIXct("2012-09-01 JST")),
            as.numeric(as.POSIXct("2014-09-30 JST"))), 
     ylim=c(0,1.2), col=rgb(0.6,0.3,0.6))
lines(sims2$time[sims2$simulation=="Stochastic"&sims2$model=="model1"], 
      sims2$a2.rev[sims2$simulation=="Stochastic"&sims2$model=="model1"],
      type="l", lwd=1.5, col=rgb(0.6,0.3,0.6))

points(dates, data$K4.NR.ACT2, lwd=1, col=rgb(0,0.7,1), pch=1, cex=0.7)
points(dates, data$K4.distNR.ACT2, lwd=1, col=rgb(0.6,0.3,0.6), pch=1, cex=0.7)

place<-c(as.POSIXct(at.dates))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]),mean(place[13:14]),mean(place[14:15]),mean(place[15:16]),mean(place[16:17]),mean(place[17:18]),mean(place[18:19]),mean(place[19:20]),mean(place[20:21]),mean(place[21:22]),mean(place[22:23]),mean(place[23:24]),mean(place[24:25]),mean(place[25:26]))
labelx<-c("9","10","11","12","1","2","3","4","5","6","7","8","9","10","11","12","1","2","3","4","5","6","7","8","9")
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
dev.off()


# Graph legends (mRNA)
pdf("figs/stochastic_legend_mRNA_200304.pdf",height=1,width=2)
par(mar=c(0,0,0,0))
par(ps=6)
par(cex=1)

plot(NA,xlim=c(0,10),ylim=c(0,10),axes=F, ann=F)

legend(1,10.5,
       legend="Obs.",lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(1,9.3,
	   legend="Sim. (stochastic, model 1, full)",
       lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,8.1,
	   legend=expression(paste("Sim. (stochastic, model 2, without ",italic(a)[N],")",sep="")),
       col="red",lty=1,pch=NA,lwd=1,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
dev.off()


# Graph legends (histone state)
pdf("figs/stochastic_legend_histone_state_200304.pdf",height=1,width=2)
par(mar=c(0,0,0,0))
par(ps=6)
par(cex=1)

plot(NA,xlim=c(0,10),ylim=c(0,10),axes=F, ann=F)

legend(1,11,
       legend=expression(paste(U[N],U[D],"(ODE)",sep="")),
       lty=1,pch=NA,lwd=0.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,9.8,
       legend=expression(paste(U[N],U[D],"(stochastic)",sep="")),
       lty=1,pch=NA,lwd=1.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,8.6,
       legend=expression(paste(M[N],U[D],"(ODE)",sep="")),
       col="red",lty=1,pch=NA,lwd=0.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,7.4,
       legend=expression(paste(M[N],U[D],"(stochastic)",sep="")),
       col="red",lty=1,pch=NA,lwd=1.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,6.2,
       legend=expression(paste(M[N],M[D],"(ODE)",sep="")),
       col="blue",lty=1,pch=NA,lwd=0.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,5,
       legend=expression(paste(M[N],M[D],"(stochastic)",sep="")),
       col="blue",lty=1,pch=NA,lwd=1.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,3.8,
       legend=expression(paste(U[N],M[D],"(ODE)",sep="")),
       col="magenta",lty=1,pch=NA,lwd=0.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,2.6,
       legend=expression(paste(U[N],M[D],"(stochastic)",sep="")),
       col="magenta",lty=1,pch=NA,lwd=1.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
dev.off()


# Graph legends (histone)
pdf("figs/stochastic_legend_histone_200305.pdf",height=1,width=2)
par(mar=c(0,0,0,0))
par(ps=6)
par(cex=1)

plot(NA,xlim=c(0,10),ylim=c(0,10),axes=F, ann=F)

legend(1,10.5,legend="Obs. (NR)",
       col=rgb(0,0.7,1),lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(1,9.3,legend="Sim. (ODE, model 1, NR)",
       col=rgb(0,0.7,1),lty=1,pch=NA,lwd=0.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,8.1,legend="Sim. (stochastic, model 1, NR)",
       col=rgb(0,0.7,1),lty=1,pch=NA,lwd=1.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,6.9,legend="Obs. (distal NR)",
       col=rgb(0.6,0.3,0.6),lty=0,pch=1,lwd=1,pt.cex=0.7,bty="n",x.intersp=0,seg.len=0.7)
legend(1,5.7,legend="Sim. (ODE, model 1, distal NR)",
       col=rgb(0.6,0.3,0.6),lty=1,pch=NA,lwd=0.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
legend(1,4.5,legend="Sim. (stochastic, model 1, distal NR)",
       col=rgb(0.6,0.3,0.6),lty=1,pch=NA,lwd=1.5,pt.cex=1,bty="n",x.intersp=0.1,seg.len=0.7)
dev.off()
