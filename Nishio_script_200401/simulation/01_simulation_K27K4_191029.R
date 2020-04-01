
################################
##### simulation_histone.R #####
################################

library('deSolve')

read.temperature <- function(csv) {
  data <- read.csv(csv)
  times <- as.POSIXct(paste(data[,'year'], data[,'month'], data[,'day'], data[,'time'], sep='/'), format='%Y/%m/%d/%H:%M')
  idx <- which(!is.na(times))
  idx <- idx[c(1, which(diff(times[idx])>0)+1)]
  idx <- idx[!is.na(data[idx, 'air.temperature'])]
  return(data.frame(time=times[idx], temperature=data[idx,'air.temperature']))
}




simulate <- function(initialState, temperature, time.start, time.step, time.end, parameters, func) {
	 parameters$func <- func
	 parameters$temperature <- approxfun(as.numeric(temperature$time - time.start), temperature$temperature)

	 times <- seq(time.start, time.end, by=time.step*60)

	 state <- ode(y = initialState, times = as.numeric(times - times[1])/60, func = modelFunc, parameters)

	temperature <- parameters$temperature(times - times[1])

	 return(list(times=times, state=state, alpha=parameters$alpha, beta=parameters$beta, gamma=parameters$gamma, epsilon=parameters$epsilon, zeta=parameters$zeta, eta=parameters$eta, iota=parameters$iota, kappa=parameters$kappa, lambda=parameters$lambda, rho=parameters$rho, phi=parameters$phi, psi=parameters$psi))

}



modelFunc <- function(t, state, parameters) {
	
  temperature <- parameters$temperature(t*60)

	return(parameters$func(state['uu'], state['mu'], state['mm'], state['um'], state['a1'], state['a2'], temperature, parameters))
}



##### models #####

model1 <- function(uu, mu, mm, um, a1, a2, temperature, parameters) {
  
  mt <- parameters$zeta / (1 + exp(parameters$alpha * (temperature - 5)))
  nt <- parameters$eta / (1 + exp(-parameters$beta * (temperature - 10)))  
  xt <- parameters$iota / (1 + exp(parameters$gamma * (temperature - 5)))
  tt <- parameters$rho / (1 + exp(-parameters$epsilon * (temperature - 15)))
  
  duu <- parameters$lambda*um - mt*uu
  dmu <- mt*uu - nt*mu
  dmm <- nt*mu - parameters$kappa*a1*mm
  dum <- parameters$kappa*a1*mm - parameters$lambda*um
  
  da1 <- tt*(1-a2)*(1-a1) - parameters$phi*(mu+mm)*a1
  da2 <- xt*(1-a1)*(1-a2) - parameters$psi*(mm+um)*a2
  
  # Convert second to day (1day=86400s)
  list(c(duu / 86400, dmu / 86400, dmm / 86400, dum / 86400, da1 / 86400, da2 / 86400))
}



model2 <- function(uu, mu, mm, um, a1, a2, temperature, parameters) {
  
  mt <- parameters$zeta / (1 + exp(parameters$alpha * (temperature - 5)))
  nt <- parameters$eta / (1 + exp(-parameters$beta * (temperature - 10)))  
  xt <- parameters$iota / (1 + exp(parameters$gamma * (temperature - 5)))
  tt <- parameters$rho / (1 + exp(-parameters$epsilon * (temperature - 15)))
  
  duu <- parameters$lambda*um - mt*uu
  dmu <- mt*uu - nt*mu
  dmm <- nt*mu - parameters$kappa*mm
  dum <- parameters$kappa*mm - parameters$lambda*um
  
  da1 <- tt*(1-a2)*(1-a1) - parameters$phi*(mu+mm)*a1
  da2 <- xt*(1-a1)*(1-a2) - parameters$psi*(mm+um)*a2
  
  # Convert second to day (1day=86400s)
  list(c(duu / 86400, dmu / 86400, dmm / 86400, dum / 86400, da1 / 86400, da2 / 86400))
}



model3 <- function(uu, mu, mm, um, a1, a2, temperature, parameters) {
  
  mt <- parameters$zeta / (1 + exp(parameters$alpha * (temperature - 5)))
  nt <- parameters$eta / (1 + exp(-parameters$beta * (temperature - 10)))  
  xt <- parameters$iota / (1 + exp(parameters$gamma * (temperature - 5)))
  tt <- parameters$rho / (1 + exp(-parameters$epsilon * (temperature - 15)))
  
  duu <- parameters$lambda*um - mt*uu
  dmu <- mt*uu - nt*mu
  dmm <- nt*mu - parameters$kappa*a1*mm
  dum <- parameters$kappa*a1*mm - parameters$lambda*um
  
  da1 <- tt*(1-a2)*(1-a1) - parameters$phi*a1
  da2 <- xt*(1-a1)*(1-a2) - parameters$psi*(mm+um)*a2
  
  # Convert second to day (1day=86400s)
  list(c(duu / 86400, dmu / 86400, dmm / 86400, dum / 86400, da1 / 86400, da2 / 86400))
}



model4 <- function(uu, mu, mm, um, a1, a2, temperature, parameters) {
  
  mt <- parameters$zeta / (1 + exp(parameters$alpha * (temperature - 5)))
  nt <- parameters$eta / (1 + exp(-parameters$beta * (temperature - 10)))  
  xt <- parameters$iota / (1 + exp(parameters$gamma * (temperature - 5)))
  tt <- parameters$rho / (1 + exp(-parameters$epsilon * (temperature - 15)))
  
  duu <- parameters$lambda*um - mt*uu
  dmu <- mt*uu - nt*mu
  dmm <- nt*mu - parameters$kappa*a1*mm
  dum <- parameters$kappa*a1*mm - parameters$lambda*um
  
  da1 <- tt*(1-a2)*(1-a1) - parameters$phi*(mu+mm)*a1
  da2 <- xt*(1-a1)*(1-a2) - parameters$psi*a2
  
  # Convert second to day (1day=86400s)
  list(c(duu / 86400, dmu / 86400, dmm / 86400, dum / 86400, da1 / 86400, da2 / 86400))
}

