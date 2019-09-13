#Dispersal Rate Calculations
  #dist: the distance dispersed in one time step (years,generations,etc.) [number]
  #time_steps: the number of time steps to model dispersal through (>0)
  #method: "mean" indicates the distance is a mean value, "median" indicates the distance is a median.
  #plot(optional): TRUE/FALSE should a plot of the total disperal probability be generated?
  #value(optional): a vector of distances.
    #If provided, will calculate the probability of dispersal to at least those distances, given time_steps

dispersalRateExp <- function(dist, time_steps, method, value, plot) {
  if (method == "mean") {
    lambda = 1 / dist
  } else if (method == "median") {
    lambda = log(2) / dist
  }

  DispFunction = function(x) {
    1 - pgamma(x, shape = time_steps, rate = lambda)
  }
  
  MeanDispersed <- time_steps / lambda
  ModeDispersed <- (time_steps - 1) / lambda
  
  if ((plot=="TRUE") | (plot == "T")) {
    par(mfrow = c(1, 2))
    XVec <- seq(from = 0, to = 2 * (time_steps / lambda), length = 1000)
    plot(XVec, dgamma(XVec, shape = time_steps, rate = lambda), 
         type = "l", 
         lwd = 2, 
         main = "Probability Density Function", 
         ylab = "Probability",
         xlab = "Distance")
    abline(v = MeanDispersed, col = "blue", lwd = 3)
    legend("topright", c(paste0("Mean = ", round(MeanDispersed, 2)), 
                        paste0("Shape = ", time_steps), 
                        paste0("Rate = ", round(lambda, 2))), 
           bty = "n",
           lwd = c(2, NA, NA), 
           col = c("blue", NA, NA))
    
    
    plot(XVec, DispFunction(XVec), 
         type = "l", 
         lwd = 2,
         main = "Cumulative Probability of Dispersal",
         ylab = "Probability",
         xlab = "Distance")
    abline(v = MeanDispersed, col = "blue", lwd = 3)
    legend("topright", c(paste0("Mean = ", round(MeanDispersed, 2))), bty = "n", lwd = c(2), col = c("blue"))
  }
    
  if (hasArg(value)){
    return(list(Value = data.frame(Distance = value, 
                                   TimeSteps = time_steps,
                                   ProbabilityofDispersal = DispFunction(value)),
                RateParameter = lambda,
                time = time_steps,
                Mean = MeanDispersed,
                Mode = ModeDispersed))
  } else {
    return(list(RateParameter = lambda, time = time_steps, Mean = MeanDispersed, Mode = ModeDispersed))
  }
}
