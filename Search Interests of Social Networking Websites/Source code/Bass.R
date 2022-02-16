# ---------------------
# Bass Diffusion Model
# ---------------------

# This file contains R functions for fitting a Bass diffusion model to data and generating forecasts from the model. The model is fit according to N(t) = m*f(t)*(1+errors), where N(t) is the demand at t, m is the number of eventual adopters, f(t) is the pdf of time to adoption in the Bass diffusion model, and (1+errors) follows a log-normal distribution. 

Bass.lik <- function(par, y, error=error) {
  # ---------------------------------------------------------------------------
  # This function calculates the log-likelihood given a set of parameters.
  # The output from this function (lik) serves as an input to the optim function.
  # Args:
  #   par --- parameters p,q,m
  #   y ---- time series
  # ----------------------------------------------------------------------------
  p <- par[1]
  q <- par[2]
  m <- par[3]
  n <- length(y)
  if (m <= 0 || p <= 0 || q <= 0 ) return(Inf) # parameter constraints
  tt <- c(1:n)
  if (error=='normal') {
    T <- c(0: n)
    F_t <- (1-exp(-(p+q)*T))/(1 + q/p*exp(-(p+q)*T))
    lik <- sum((m*(F_t[2: length(T)] - F_t[1: (length(T) - 1)]) - y)^2)
  } else {
    lik <- sum((log(m*(p+q)^2/p*(exp(-(p+q)*tt))/(1 + q/p*exp(-(p+q)*tt))^2) - log(y))^2)
  }
  return(lik)
}

Bass.fit <- function(y, method = "Nelder-Mead", error="lognormal") {
  # ------------------------------------------------------------------
  # This function fits the Bass Diffusion Model using MLE.
  # Output:
  #   model parameters,in-sample fitted values,peak time and peak value.
  # Args:
  #   y ---- time series
  #   method ---- method argument passed to optim
  #   error ---- error type; can be set to "normal" or "lognormal"
  # -------------------------------------------------------------------
  # Set initial parameters for MLE
  p <- 0.01; q <- 0.1; m <- sum(y)
  initial_value <- c(p,q,m)
  # Estimate parameters with MLE using the optim function
  model <- optim(par = initial_value, fn = Bass.lik, y = y,
                 method = "Nelder-Mead", error=error,
                 control = list(maxit = 2000))
  parameter <- model$par
  model <- optim(par = parameter, fn = Bass.lik, Bass.lik, y = y,
                 method = "Nelder-Mead",error=error,
                 control = list(maxit = 2000 ))
  parameter  <- model$par
  
  # Calculate in-sample fitted values
  n <- length(y)
  p <- parameter[1]
  q <- parameter[2]
  m <- parameter[3]
  tt <- c(1:length(y))
  if (error=='normal') {
    T <- c(0:n)
    fit.F <- (1-exp(-(p+q)*T))/(1 + q/p*exp(-(p+q)*T))
    fit.y <- m*(fit.F[2: length(T)] - fit.F[1: (length(T) - 1)])
  } else {
    fit.y <- m*(p+q)^2/p*(exp(-(p+q)*tt))/(1 + q/p*exp(-(p+q)*tt))^2 
  }

  sumy <- cumsum(y)
  if (error=='normal') {
    sigma <- sqrt(sum((fit.y-y)^2)/length(y))
  } else {
    sigma <- sqrt(sum((log(fit.y)-log(y))^2)/length(y)) # calculate sd of error terms
  }

  # Calculate peak time and peak value
  mode <- -log(p/q)/(p+q) # calculate peak time
  if (mode < 1) mode <- ceiling(-log(p/q)/(p+q))
  if (q > p) {
    peakValue <- m*(p+q)^2/(4*q) # calculate peak value
  } else {
    peakValue <- m*p
  }
  output <- list(par = parameter, fit.y = fit.y, y = y, 
                 peakValue = peakValue, peakTime = mode, sigma = sigma)
  return(structure(output, class = "Bass.model"))
}

forecast.Bass <- function(model,h,PI = TRUE, level=c(80, 95), error='lognormal') {
  # -------------------------------------------------------------------
  # This function generates forecasts from the Bass model.
  # Args:
  #   model ---- the Bass Diffusion model
  #   h ---- forecasting horizon
  #   PI ---- if FALSE only point forecasts will be generated
  #   level --- level of prediction intervals;default is 80% and 95%.
  #   error ---- error type; can be set to "normal" or "lognormal"
  # --------------------------------------------------------------------
  sigma <- model$sigma
  parameter <- model$par  
  p <- parameter[1]
  q <- parameter[2]
  m <- parameter[3]
  n <- length(model$y)
  if (error=='normal') {
    T <- (n ) : (n + h)
    fit.F <- (1-exp(-(p+q)*T))/(1 + q/p*exp(-(p+q)*T))
    yf <- m*(fit.F[2: length(T)] - fit.F[1: (length(T) - 1)])
  } else {
    tt <- (n+1) : (n+h)
    yf <- m*(p+q)^2/p*(exp(-(p+q)*tt))/(1 + q/p*exp(-(p+q)*tt))^2
  }
  y <- model$y
  
  if (PI==TRUE) {
    interval <- matrix(0, h, 2*length(level))
    name <- matrix(0, 1, 2*length(level))
    for ( i in 1:length(level) ) {
      name[2*i] <- paste("Hi", level[i])
      name[2*(i-1)+1] <- paste("Lo", level[i])
    }
    # calculate quantiles according to the log-normal distribution.
    for (i in 1:length(level)){ 
      if (error=='normal') {
        interval[, 2*(i-1)+1] <-  qnorm(0.5 - level[i]/200, mean = yf, sd = sigma)
        interval[, 2*i] <- qnorm(0.5 + level[i]/200, mean = yf, sd = sigma)
      } else {
        interval[, 2*(i-1)+1] <-  qlnorm(0.5 - level[i]/200, mean = log(yf), sd = sigma)
        interval[, 2*i] <- qlnorm(0.5 + level[i]/200, mean = log(yf), sd = sigma) 
      }
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

plot.forecast.Bass <- function(prediction,PI = TRUE,ylim = NULL,main=NULL,xlab=NULL,ylab=NULL){
  # ------------------------------------------------------------------------
  # This function creates a plot of the forecasts generated from a Bass model.
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
      title <- 'Forecasts from Bass'
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
    plot(y, xlim = c(1,(length(y)+h)), ylim = ylim, type = 'l',
         main = title,xlab = xlab,ylab = ylab,axes = FALSE)
    
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

