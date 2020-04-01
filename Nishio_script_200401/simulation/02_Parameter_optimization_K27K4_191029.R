
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

# Index for sampling time
data.idx <- as.numeric((dates - time.start) / time.step, units='mins')

# Trimming of temperature data
temperature <- temperature[time.start<=temperature$time & 
				temperature$time<=time.end+as.difftime(1,units='days'), ]

# Setting of parameter grid
grid.parameters <- as.data.frame(
  list(
     alpha=runif(1000,0,100),
     beta=runif(1000,0,100),
     gamma=runif(1000,0,100),
	 epsilon=runif(1000,0,100),
	 zeta=runif(1000,0,10),
     eta=runif(1000,0,10),
     iota=runif(1000,0,10),
     kappa=runif(1000,0,10),
	 lambda=runif(1000,0,10),
	 rho=runif(1000,0,10),
	 phi=runif(1000,0,10),
	 psi=runif(1000,0,10)
  )
)

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

# Comparison between simulation and observation
evaluate.param <- function(
  model.fun, 
  parameter, 
  temperature, 
  data, data.idx,
  time.start, time.step, time.end, 
  init=NULL){
  
  sim.results <- try(simulate.param(
    model.fun=model.fun,
    parameter=parameter,
    temperature=temperature,
    time.start=time.start,
    time.step=time.step,
    time.end=time.end,
    init=init
    ))
      # Evaluate the results if simulation is successful
      if (
        class(sim.results) != "try-error" &&
        sum(is.na(sim.results$state))==0 && 
      	dim(sim.results$state)[1] == length(sim.results$times) && 
      	dim(sim.results$state)[2] == 7){
        
        if((sim.results$alpha > 100) || (sim.results$beta > 100) ||
        (sim.results$gamma > 100) || (sim.results$epsilon > 100) ||
        (sim.results$zeta > 100) || (sim.results$eta > 100) || 
        (sim.results$iota > 100) || (sim.results$kappa > 100) ||
        (sim.results$lambda > 100) || (sim.results$rho > 100) || 
        (sim.results$phi > 100) || (sim.results$psi > 100)){
          value <- Inf
        }
        else{
        
          value <-
          	sum(((max(data$K27.NR.STM)*(sim.results$state[data.idx, 'mu'] + 
          		sim.results$state[data.idx, 'mm']) -
          		data$K27.NR.STM)/sd(data$K27.NR.STM))^2)
          value <- value + 
          	sum(((max(data$K27.distNR.STM)*(sim.results$state[data.idx, 'um'] + 
          		sim.results$state[data.idx, 'mm']) - 
          		data$K27.distNR.STM)/sd(data$K27.distNR.STM))^2)
          value <- value + 
          	sum(((max(data$K4.NR.ACT2)*sim.results$state[data.idx, 'a1'] - 
          		data$K4.NR.ACT2)/sd(data$K4.NR.ACT2))^2)   
          value <- value + 
          	sum(((max(data$K4.distNR.ACT2)*sim.results$state[data.idx, 'a2'] - 
          		data$K4.distNR.ACT2)/sd(data$K4.distNR.ACT2))^2)   
      	}
      }
      else{
          value <- Inf
  }
  return(value)
}

save.image("Model_191029.RData")





# model1
setwd("~/Documents/nishio/simulation/")
source('01_simulation_K27K4_191029.R')
load('Model_191029.RData')

# Evaluation of each row of the parameter grid (SANN)
calc.time1<-system.time(
grid.cor <- apply(grid.parameters, 1, function(g)
  evaluate.param(
    model1,
    as.list(g),
    temperature,
    data,
    data.idx,
    time.start,
    time.step,
    time.end
    )
  )
)

# Further optimization using SANN method
grid.best <- grid.parameters[which.min(grid.cor),]
print(grid.best)

calc.time2<-system.time(
optim.result <- optim(
  log(grid.best), 
  function(p) 
    evaluate.param(
      model1,
      as.list(exp(p)),
      temperature,
      data,
      data.idx,
      time.start,
      time.step,
      time.end
      ), 
  method="SANN",
  control=list(trace=1, maxit=5000)
  )
)

parameter.optim <- exp(optim.result$par)
parameter.optim

script <- 'parameter_optimization_K27K4.R'
save(parameter.optim, grid.best, grid.parameters, grid.cor, script, 
		file='optim/optim_parameters_K27K4_model1_SANN_191029')




# model2
# Evaluation of each row of the parameter grid (SANN)

setwd("~/Documents/nishio/simulation/")
source('01_simulation_K27K4_191029.R')
load('Model_191029.RData')

calc.time1<-system.time(
grid.cor <- apply(grid.parameters, 1, function(g)
  evaluate.param(
    model2,
    as.list(g),
    temperature,
    data,
    data.idx,
    time.start,
    time.step,
    time.end
    )
  )
)

# Further optimization using SANN method
grid.best <- grid.parameters[which.min(grid.cor),]
print(grid.best)

calc.time2<-system.time(
optim.result <- optim(
  log(grid.best), 
  function(p) 
    evaluate.param(
      model2,
      as.list(exp(p)),
      temperature,
      data,
      data.idx,
      time.start,
      time.step,
      time.end
      ), 
  method="SANN",
  control=list(trace=1, maxit=5000)
  )
)

parameter.optim <- exp(optim.result$par)
parameter.optim

script <- 'parameter_optimization_K27K4.R'
save(parameter.optim, grid.best, grid.parameters, grid.cor, script, 
		file='optim/optim_parameters_K27K4_model2_SANN_191029')




# model3
# Evaluation of each row of the parameter grid (SANN)

setwd("~/Documents/nishio/simulation/")
source('01_simulation_K27K4_191029.R')
load('Model_191029.RData')

calc.time1<-system.time(
grid.cor <- apply(grid.parameters, 1, function(g)
  evaluate.param(
    model3,
    as.list(g),
    temperature,
    data,
    data.idx,
    time.start,
    time.step,
    time.end
    )
  )
)

# Further optimization using SANN method
grid.best <- grid.parameters[which.min(grid.cor),]
print(grid.best)

calc.time2<-system.time(
optim.result <- optim(
  log(grid.best), 
  function(p) 
    evaluate.param(
      model3,
      as.list(exp(p)),
      temperature,
      data,
      data.idx,
      time.start,
      time.step,
      time.end
      ), 
  method="SANN",
  control=list(trace=1, maxit=5000)
  )
)

parameter.optim <- exp(optim.result$par)
parameter.optim

script <- 'parameter_optimization_K27K4.R'
save(parameter.optim, grid.best, grid.parameters, grid.cor, script, 
		file='optim/optim_parameters_K27K4_model3_SANN_191029')




# model4
# Evaluation of each row of the parameter grid (SANN)

setwd("~/Documents/nishio/simulation/")
source('01_simulation_K27K4_191029.R')
load('Model_191029.RData')

calc.time1<-system.time(
grid.cor <- apply(grid.parameters, 1, function(g)
  evaluate.param(
    model4,
    as.list(g),
    temperature,
    data,
    data.idx,
    time.start,
    time.step,
    time.end
    )
  )
)

# Further optimization using SANN method
grid.best <- grid.parameters[which.min(grid.cor),]
print(grid.best)

calc.time2<-system.time(
optim.result <- optim(
  log(grid.best), 
  function(p) 
    evaluate.param(
      model4,
      as.list(exp(p)),
      temperature,
      data,
      data.idx,
      time.start,
      time.step,
      time.end
      ), 
  method="SANN",
  control=list(trace=1, maxit=5000)
  )
)

parameter.optim <- exp(optim.result$par)
parameter.optim

script <- 'parameter_optimization_K27K4.R'
save(parameter.optim, grid.best, grid.parameters, grid.cor, script, 
		file='optim/optim_parameters_K27K4_model4_SANN_191029')

