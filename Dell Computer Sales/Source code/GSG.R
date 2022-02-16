# -----------------------------
# Gamma/shifted Gompertz Model
# -----------------------------

# This file contains R functions for fitting a Gamma/shifted Gompertz model to data and generating forecasts from the model. The model is fit according to N(t) = m*f(t)*(1+errors), where N(t) is the demand at t, m is the number of eventual adopters, f(t) is the pdf of time to adoption in the Gama/shifted Gompertz model, and (1+errors) follows a log-normal distribution.

GSG.lik <- function(par, y){
  # ---------------------------------------------------------------------------
  # This function calculates the likelihood function given a set of parameters.
  # The output from this function (lik) serves as an input to the optim function.
  # Args:
  #   par --- parameters m,lambda,nu,mu
  #   y ---- time series
  # ----------------------------------------------------------------------------
  m <- par[1]
  lambda <- par[2]
  nu <- par[3]
  mu <- par[4]
  n <- length(y)
  if(m < 0 | (nu) < 0 | mu < 0 | lambda < 0 ) return(Inf) # parameter constraints
  tt <- c(1:length(y))
  lik <- sum((log(m*lambda/nu*(exp(-lambda*tt)*(mu+nu+(1-mu)*exp(-lambda*tt)))/(1+1/nu*exp(-lambda*tt))^(mu+1)) - log(y))^2)
  return(lik)
} 

GSG.fit <- function(y=NULL, method = "Nelder-Mead"){
  # ------------------------------------------------------------------
  # This function fits the Gamma/shifted Gompertz Model using MLE.
  # Output:
  #   model parameters,in-sample fitted values,peak time and peak value.
  # Args:
  #   y ---- time series
  #   method ---- method argument passed to optim
  # -------------------------------------------------------------------
  # Set initial parameters for MLE
  m <- sum(y); lambda <- 0.11; nu <- 0.1; mu <- 1
  initial <- c(m, lambda, nu, mu)
  # Estimate parameters with MLE using the optim function
  modelSG <- optim(par = initial, fn = GSG.lik, y = y, 
                   method = method, control = list(maxit = 2000))
  parameter <- modelSG$par
  modelSG <- optim(par = parameter, fn = GSG.lik, y = y, 
                   method = method, control = list(maxit = 2000))
  
  # Calculate in-sample fitted values
  m <-  parameter[1]
  lambda <- parameter[2]
  nu <- parameter[3]
  mu <- parameter[4]
  n <- length(y)
  tt <- c(1:length(y))
  fit.y <- m*lambda/nu*(exp(-lambda*tt)*(mu+nu+(1-mu)*exp(-lambda*tt)))/(1+1/nu*exp(-lambda*tt))^(mu+1)
  sigma <- sqrt(sum((log(fit.y)-log(y))^2)/length(y))
  
  # Calculate peak time and peak value
  A <- -(1/nu)*(mu-1)^2
  B <- (1/nu)*mu^2 + 3*mu - 2
  C <- - 1/ (1/nu) - mu
  D <- mu*(-4 + 5*mu - 4*(1/nu) + 4*mu*(1/nu) + 2*mu^2*(1/nu) + mu^3*(1/nu)^2)
  if (B < 0 || D < 0) {
    peakValue <- lambda*(nu/(1+nu))^mu *m
  } else {
    root <- 2*C/(-B-sqrt(B^2 - 4*A*C))
    root_2 <- (-B - sqrt(B^2 - 4*A*C))/(2*A) 
    if (mu == 1) {
      mode <- -log(nu)/lambda
      if (mode < 1) mode <- ceiling(-log(root)/lambda)
      mode1 <- mode - 1
      if (mode <= 0) {
        mode = 0
        mode1 = 1
      }
      F_mode <- (1-exp(-lambda*mode))/((1+1/nu*exp(-lambda*mode))^mu)
      F_mode1 <- (1-exp(-lambda*mode1))/((1+1/nu*exp(-lambda*mode1))^mu)  
      peakValue <- m*(F_mode - F_mode1)
    } else {
      if (root > 0 & root < 1 & root_2 >= 1){
        mode <- -log(root)/lambda
        if (mode < 1) mode <- ceiling(-log(root)/lambda)
        mode1 <- mode - 1
        if (mode < 1) mode <- ceiling(-log(root)/lambda)
        mode1 <- mode - 1
        if (mode <= 0) {
          mode = 0
          mode1 = 1
        }
        F_mode <- (1-exp(-lambda*mode))/((1+1/nu*exp(-lambda*mode))^mu)
        F_mode1 <- (1-exp(-lambda*mode1))/((1+1/nu*exp(-lambda*mode1))^mu)  
        peakValue <- m*(F_mode - F_mode1)
      } else {
        if (root > 0 & root < 1 & root_2 < 1 & root_2 > 0) {
          mode0 <- 0
          mode <- -log(root)/lambda
          if (mode < 1) mode <- ceiling(-log(root)/lambda)
          mode1 <- mode - 1
          if (mode < 1) mode <- ceiling(-log(root)/lambda)
          mode1 <- mode - 1
          if (mode <= 0) {
            mode = 0
            mode1 = 1
          }
          peak0 <- lambda*(nu/(1+nu))^mu
          F_mode <- (1-exp(-lambda*mode))/((1+1/nu*exp(-lambda*mode))^mu)
          F_mode1 <- (1-exp(-lambda*mode1))/((1+1/nu*exp(-lambda*mode1))^mu)  
          peakValue <- m*(F_mode - F_mode1)
          peakValue <- max(peakValue, peak0)
        } else {
          peakValue <- lambda*(nu/(1+nu))^mu
        }
      }
    }
  }
  out <- list(par = parameter, fit.y = fit.y, y = y, 
              sigma = sigma, peakTime = mode, peakValue = peakValue)
  return(structure(out, class = "GSG"))
}

forecast.GSG <- function(model,h, PI = TRUE, level=c(80, 95)) {
  # -------------------------------------------------------------------
  # This function generates forecasts from a GSG model.
  # Args:
  #   model ---- the GSG model
  #   h ---- forecasting horizon
  #   PI ---- if FALSE only point forecasts will be generated
  #   level --- level of prediction intervals;default is 80% and 95%.
  # --------------------------------------------------------------------
  sigma <- model$sigma
  par <- model$par  
  m <- par[1]
  lambda <- par[2]
  nu <- par[3]
  mu <- par[4]
  n <- length(model$y)

  tt <- c(n+1):c(n+h)
  yf <- m*lambda/nu*(exp(-lambda*tt)*(mu+nu+(1-mu)*exp(-lambda*tt)))/(1+1/nu*exp(-lambda*tt))^(mu+1)
  y <- model$y
  
  if (PI==TRUE) {
    interval <- matrix(0, h, 2*length(level))
    name <- matrix(0, 1, 2*length(level))
    for ( i in 1:length(level) ) {
      name[2*i] <- paste("Hi", level[i])
      name[2*(i-1)+1] <- paste("Lo", level[i])
    }
    for (i in 1:length(level)){
      interval[, 2*(i-1)+1] <-  qlnorm(0.5 - level[i]/200, mean = log(yf), sd = sigma)
      interval[, 2*i] <- qlnorm(0.5 + level[i]/200, mean = log(yf), sd = sigma)
    }
    colnames(interval) <- name
  }
  if (PI == TRUE){
    out <- cbind(yf, interval)
  } else {
    out <- yf
  }
  result <- list(y = model$y, 'predictions' = out, 
                 edian = yf, level=level, PI = PI)
  return(structure(result, class = "GSG.model")) 
}
