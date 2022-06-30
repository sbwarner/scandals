# Functions to plot acf and pacf with ggplot2


#  a function to plot the ACF using ggplot2
ggacf <- function(y, lag = 12, plot.zero="no", alpha = 0.05)
{
  #y<- y.dt
  #plot.zero="no"
  #alpha=0.05
  T    <- length(y)
  y.acf <- acf(y, lag.max=lag, plot=FALSE)
  if (plot.zero == "yes") y.acf <- data.frame(lag= 0:lag, acf = y.acf$acf)
  if (plot.zero == "no")  y.acf <- data.frame(lag= 1:lag, acf = y.acf$acf[-1])
  
  library(ggplot2)
  ggplot(y.acf, aes(lag, acf)) + 
    geom_bar(stat="identity", fill="orange") + 
    geom_hline(yintercept =  qnorm(1-alpha/2) * T^(-0.5), color="steelblue3", linetype="dashed") + 
    geom_hline(yintercept =  qnorm(alpha/2) * T^(-0.5), color="steelblue3", linetype="dashed") +
    geom_hline(yintercept = 0, color="steelblue3") +   
    theme_classic() + ylab("") + ggtitle("ACF")
}
ggpacf <- function(y, lag = 12, plot.zero="no", alpha = 0.05)
{
  T    <- length(y)

  #plot.zero="yes"
   y.pacf <- pacf(y, lag.max=lag, plot=FALSE)
   if (plot.zero == "yes")  {
    lag0pacf <- 1
  toplotpacf <- c(1, y.pacf$acf[-1])
  toplotlags <- 0:length(y.pacf$lag[-1])
  y.pacf <- data.frame(lag= toplotlags, pacf = toplotpacf)}
  if (plot.zero == "no") y.pacf <- data.frame(lag= 1:lag, pacf = y.pacf$acf)
  
  
  
  
  library(ggplot2)
  ggplot(y.pacf, aes(lag, pacf)) + 
    geom_bar(stat="identity", fill="orange") + 
    geom_hline(yintercept =  qnorm(1-alpha/2) * T^(-0.5), color="steelblue3", linetype="dashed") + 
    geom_hline(yintercept =  qnorm(alpha/2) * T^(-0.5), color="steelblue3", linetype="dashed") +
    geom_hline(yintercept = 0, color="steelblue3") +   
    theme_classic() + ylab("") + ggtitle("PACF")
}
