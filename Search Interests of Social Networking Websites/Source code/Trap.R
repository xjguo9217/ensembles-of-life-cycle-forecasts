# -----------------
# Trapezoid Model
# -----------------

# This file contains R functions for fitting a Trapezoid model to data and generating forecasts from the model. The model is fit according to N(t) = m*f(t)*(1+errors), where N(t) is the demand at t, m is the number of eventual adopters, f(t) is the pdf of time to adoption in the Trapezoid model, and (1+errors) follows a log-normal distribution. 

Trapezoid.lik <- function(par, y) {
  # ---------------------------------------------------------------------------
  # This function calculates the likelihood function given a set of parameters.
  # The output from this function (lik) serves as an input to the optim function.
  # Args:
  #   par --- parameters p,q,m
  #   y ---- time series
  # ----------------------------------------------------------------------------
  a <- par[1]
  b <- par[2]
  c <- par[3]
  tau1 <- par[4]
  deltatau <- par[5]
  tau2 <- tau1 + deltatau
  if (tau1 > tau2|c > 0|deltatau < 0|a < 0) return(Inf)
  fit.y <- rep(NA, length(y))
  for (t in 1: length(y)) {
    if (t < tau1) fit.y[t] <- a*t + b
    if (t >= tau1 & t < tau2) fit.y[t] <- a*tau1 + b
    if (t >= tau2) fit.y[t] <- c*(t - tau2) + (a*tau1 + b)
  }
  fit.y[fit.y <= 0] <- min(fit.y[fit.y>0])/10
 lik <- sum((log(fit.y) - log(y))^2)
  return(lik)
}

Trapezoid.fit <- function(y, method = "Nelder-Mead") {
  # ------------------------------------------------------------------
  # This function fits the Trapezoid Model using MLE.
  # Output:
  #   model parameters and in-sample fitted values.
  # Args:
  #   y ---- time series
  #   method ---- method argument passed to optim
  # -------------------------------------------------------------------
  # Set initial parameters for MLE
  ylm <- y[1:which.max(y)]
  if (length(ylm) < 5) ylm <- y[1:5]
  aaa=lsfit(1:length(ylm),ylm)
  ini_a <- aaa$coefficients[2]
  ini_b <- aaa$coefficients[1]
  ini_tau1 <- which.max(y)
  deltatau <- 0
  ylm <- y[which.max(y):length(y)]
  if (length(ylm) < 5) ylm <- y[(length(y)-5):length(y)]
  aaa=lsfit(1:length(ylm),ylm)
  ini_c <- aaa$coefficients[2]
  if (ini_a < 0) ini_a <- 0.01
  if (ini_c > 0) ini_c <- -ini_a

  initial_value <- c(ini_a, ini_b, ini_c, ini_tau1, deltatau)
  model <- optim(par = initial_value, fn = Trapezoid.lik, y = y,
                 method = "Nelder-Mead",
                 control = list(maxit = 2000))
  parameter <- model$par
  model <- optim(par = parameter, fn = Trapezoid.lik, y = y,
                 method = "Nelder-Mead",
                 control = list(maxit = 2000))
  
  # Calculate in-sample fitted values
  par <- model$par
  a <- par[1]
  b <- par[2]
  c <- par[3]
  tau1 <- par[4]
  deltatau <- par[5]
  tau2 <- tau1 + deltatau
  
  fit.y <- rep(NA, length(y))
  for (t in 1: length(y)) {
    if (t < tau1) fit.y[t] <- a*t + b
    if (t >= tau1 & t < tau2) fit.y[t] <- a*tau1 + b
    if (t >= tau2) fit.y[t] <- c*(t - tau2) + (a*tau1 + b)
  }
  fit.y[fit.y <= 0] <- min(fit.y[fit.y>0])/10
  sigma <- sqrt(mean(((log(y) - log(fit.y)))^2))
  output <- list(par = par, fit.y = fit.y, y = y, sigma = sigma )
  return(structure(output, class = "Trapezoid.model"))
}

forecast.Trapezoid <- function(model, h, PI = TRUE, level=c(80, 95)) {
  # -------------------------------------------------------------------
  # This function generates forecasts from a Trapezoid model.
  # Args:
  #   model ---- the Trapezoid model
  #   h ---- forecasting horizon
  #   PI ---- if FALSE only point forecasts will be generated
  #   level --- level of prediction intervals;default is 80% and 95%.
  # --------------------------------------------------------------------
  sigma <- model$sigma
  par <- model$par  
  a <- par[1]
  b <- par[2]
  c <- par[3]
  tau1 <- par[4]
  tau2 <- par[4] + par[5]
  n <- length(model$y)
  T <- (n) : (n + h)
  y <- model$y
  
  yf <- rep(NA, h)
  for (t in T) {
    if (t < tau1) yf[t-n] <- a*t + b
    if (t >= tau1 & t < tau2) yf[t-n] <- a*tau1 + b
    if (t >= tau2) yf[t-n] <- c*(t - tau2) + (a*tau1 + b)
  }
  n_0 <- which(yf <= 0)
  yf[yf <= 0] <- min(yf[yf>0])/10
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
                 median = yf, level=level, PI = PI)
  return(structure(result, class = "Bass.model")) 
}

plot.forecast.Trapezoid <- function(prediction,PI = TRUE,ylim = NULL,main=NULL,xlab=NULL,ylab=NULL){
  # ------------------------------------------------------------------------
  # This function creates a plot of the forecasts generated from a Trap model.
  # Args:
  #   prediction ---- prediction form the model
  #   PI ---- if FALSE only point forecasts will be plotted
  #   ylim, main, xlab, ylab --- graphic parameters
  # -------------------------------------------------------------------------
  y <- prediction$y
  forecasts <- prediction$mean
  if (PI == FALSE) {
    if (is.null(xlim))
      xlim <- c(0, (length(y) + length(forecasts)))
    if (is.null(ylim))
      ylim <- c(min(y, forecasts), max(y, forecasts))
    plot(y, main = main, xlim = xlim, ylim = ylim, type = 'l')
    if (length(forecasts) == 1) {
      points((length(y)+1), forecasts, col = 'blue', lwd = 2, pch=19)
    } else {
      lines((length(y)+1):(length(y)+length(forecasts)), forecasts, col = 'blue')
    }
    lines(pointForecast$fit, col = 'red')
  } else {
    level <- prediction$level
    h <- length(prediction$prediction[,1])
    out <- prediction$predictions
    if(any(is.na(out)))
      PI <- FALSE
    if (h == 1){
      out <- out[-c(2,3)]
    }else{
      if(PI == TRUE){
        out <- cbind(out[,1],out[,2:ncol(out)])
      } else {
        out <- matrix(out[,1])
      }
    }
    if (is.null(main)){
      title <- 'Forecasts from Trapezoid'
    }else{
      title <- main
    }
    if (is.null(ylim) ) {
      ylim <- c(min(min(out[!is.na(out)]), y), 
                max(max(out[!is.na(out)]),y))
      if (ylim[2] == Inf){
        n0 <- which(out == Inf)
        ylim[2] <- max(out[-n0])
      }
    }
    plot(y, xlim = c(1,(length(y)+h)), ylim = ylim, type = 'l',main = title,xlab = xlab,ylab = ylab,axes = FALSE)
    
    xxx <- (length(y)+1):(length(y)+h)
    if (h == 1){
      points(xxx,out[1],col = 'blue',lwd = 2)
    } else {
      lines(xxx,out[,1],col = 'blue',lwd = 2,xlab = xlab)
    }
    
    for (i in 1:length(level)){
      if (h == 1){
        out1 <- out[2:length(out)]
      }else{
        out1 <- out[,2:(dim(out)[2])]
      }
      j <- rev(1:(2*length(level)))
      shadecols <- rev(colorspace::sequential_hcl(52 )[31]) #[level - 49]
      if (h == 1){
        polygon(c(xxx + c(-0.5,0.5,0.5,-0.5)), 
                c(rep(out1[c(j[(i-1)*2+1])],2),rep(out1[c(j[i*2])], 2)), 
                col  =  shadecols[i], border  =  FALSE) 
        points(xxx,out[1],col = 'blue',lwd = 2, pch=19)
      } else {
        polygon(c(xxx, rev(xxx)), 
                c(out1[, c(j[(i-1)*2+1])],rev(out1[, c(j[i*2])])), 
                col  =  shadecols[i], border  =  FALSE)
        lines(xxx,out[,1],col = 'blue',lwd = 2,xlab = xlab)
      }
    }
  }
  axis(1, pos=0) 
  axis(2, pos=0) 
}


