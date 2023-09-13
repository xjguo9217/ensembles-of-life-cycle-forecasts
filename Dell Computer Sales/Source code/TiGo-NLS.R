# -------------------------------
# Tilted-Gompertz Diffusion Model
# -------------------------------

# This file contains R functions for fitting a tilted-Gompertz diffusion model to data and generating forecasts from the model. The model is fit according to the TiGo formulation rather than the exponential smoothing formulation, which is N(t) = m*f(t)*(1+errors), where N(t) is the demand at t, m is the number of eventual adopters, f(t) is the pdf of time to adoption in the Bass diffusion model, and (1+errors) follows a log-normal distribution. 

# The gamma function in R base will return Inf when x in gamma(x) is greater than 171. Therefore, we need the Rmpfr function to calculate the gamma function and lower incomplete gamma function when x > 171.
library(Rmpfr)
gammaMpfr <- function(x){
  if (x > 170){
    return( gamma(as(x,"mpfr")))
  }else{
    return(gamma(x))
  }
}

lower.inc.gamma <- function(a, x) {
  lowerIncGamma <- pgamma(x, a, lower = TRUE) * gammaMpfr(a)[[1]]
  return(lowerIncGamma)
}

TiGo.lik <- function(par, y, skew, maxPeak, ultimateAdopt,lambda, rho, delta, m, error=error){
  # ---------------------------------------------------------------------------
  # This function calculates the log-likelihood given a set of parameters.
  # The output from this function (lik) serves as an input to the optim function.
  # Args:
  #   par --- parameters m,lambda,delta,rho
  #   y ---- time series
  # ----------------------------------------------------------------------------
  fixed <- c(NA, NA, NA, NA)
  coef1 <- c(0,0,0,0,0)
  if (!is.null(m))
    fixed[1] <- m
  if (!is.null(lambda))
    fixed[2] <- lambda
  if (!is.null(delta))
    fixed[3] <- delta
  if (!is.null(rho))
    fixed[4] <- rho
  if ( all(is.na(fixed)) == TRUE )
    fixed=NULL
  if (is.null(fixed) == FALSE){
    nn  <- which(is.na(fixed))
    nn1 <- which(!is.na(fixed))    
    coef1[nn1] <- fixed[nn1]
    coef1[nn]  <- par
  }else{
    coef1 <- par  
  }
  
  m <- coef1[1]
  lambda <- coef1[2]
  delta <- coef1[3]
  rho <- coef1[4]
  if (!is.null(ultimateAdopt) && m > ultimateAdopt) return(Inf)
  if (skew == 'R' & lambda <= 0 ) return(Inf) 
  if (skew == 'L' & lambda >= 0 ) return(Inf) 
  # parameter constraints
  if (rho <= 0 || m <= 0 || delta <= 0) return(Inf) 
  n <- length(y)
  T <- c(1): (n)
  I <- ifelse(lambda < 0, 1, 0)
  if (delta > 171.6) return(Inf) 
  lowerGamma <- pgamma(rho, delta, lower = TRUE) * gamma(delta)
  gammaF <- gamma(delta) 
  if (lambda*(rho - delta) > 0) {
    log_peak <- log(m) + log(abs(lambda))- log(abs(lowerGamma - I * gammaF)) + delta * log(delta) - delta
  } else {
    log_peak <- log(m) - delta + delta*log(delta)+ log(abs(lambda))  - log(abs(lowerGamma - I * gammaF))  
  }
  if (is.null(maxPeak)) {
    if (is.na(log_peak)) return(Inf)
  } else {
    if (is.na(log_peak) || as.numeric(log_peak) > log(maxPeak)) return(Inf)
  }
  yf <- m * abs(lambda) * exp(-lambda*T*delta + delta*log(rho)-rho*exp(-lambda*T)-log(abs(lowerGamma - I * gammaF)))
  if (error == 'normal') {
    sse <- sum((y-yf)^2)
  } else {
    sse <- sum((log(y)-log(yf))^2)
  }
  return(as.numeric(sse))
} 

TiGo.NLS <- function(y=NULL, maxPeak = NULL, ultimateAdopt = NULL,lambda=NULL, rho=NULL, delta=NULL, M = NULL, control = list(maxit = 20000), method = "Nelder-Mead", error="lognormal"){
  # ---------------------------------------------------------------------------
  # This function fits the TiGo Diffusion Model using MLE.
  # Output:
  #   model parameters,in-sample fitted values,peak time and peak value.
  # Args:
  #   y ---- time series
  #   maxPeak ---- fix a peak value if known
  #   ultimateAdopt ---- fix the ultimate adoption (market potential) if known
  #   method ---- method argument passed to optim
  #   error ---- error type; can be set to "normal" or "lognormal"
  # ----------------------------------------------------------------------------
  fixed <- c(NA, NA, NA, NA)
  if (!is.null(M))
    fixed[1] <- M
  if (!is.null(lambda))
    fixed[2] <- lambda
  if (!is.null(delta))
    fixed[3] <- delta
  if (!is.null(rho))
    fixed[4] <- rho
  mi <- M
  # Fit the right skewed model
  # initialize parameters
  if (is.null(M)) {
    mi <- sum(y)
    if (!is.null(ultimateAdopt) && mi>ultimateAdopt) {
      mi <- ultimateAdopt - 1e-5
    }
  }
  
  if (is.null(lambda)) {
    lambdai <- 0.03
  }else{
    lambdai <- lambda
  }
  if (is.null(delta)) {
    deltai <- 0.67
  } else {
    deltai <- delta
  }
  if (is.null(rho)) {
    rhoi <- 3.74
  }else {
    rhoi <- rho
  }
 
  c <- lambdai*rhoi^deltai/(lower.inc.gamma(deltai,rhoi)[1])
  if (!is.null(maxPeak) && mi*c*(deltai/rhoi)^deltai*exp(-deltai) > maxPeak) {
    lambdai <- maxPeak/mi*exp(deltai)*deltai^(-deltai)*lower.inc.gamma(deltai,rhoi)[1]
  }
  initial.valuer <- c(sum(y), lambdai, deltai, rhoi)
  set <- which(is.na(fixed))
  initial.valuer <- initial.valuer[set]
    modelTG <- optim(par = initial.valuer, fn = TiGo.lik, skew = 'R', maxPeak = maxPeak, error=error,
                     y = y, ultimateAdopt = ultimateAdopt, m = M, lambda=lambda, delta = delta, rho = rho,
                     method = method, control = control)
    initial <- modelTG$par
    modelTG <- optim(par = initial, fn = TiGo.lik, skew = 'R', maxPeak = maxPeak, error=error,
                     y = y, ultimateAdopt = ultimateAdopt,m = M, lambda=lambda, delta = delta, rho = rho,
                     method = method, control = control)
  
  # Fit the left skewed model
  if (is.null(lambda)) lambdai <- -0.03
  if (is.null(delta)) deltai <- 0.67
  if (is.null(rho)) rhoi <- 0.17
  
  lowerGamma <- pgamma(rhoi, deltai, lower = TRUE) * gamma(deltai)
  gammaF <- gamma(deltai) 
  c <- lambdai*rhoi^deltai/(lowerGamma - gammaF )
  if (!is.null(maxPeak) && mi*c*(deltai/rhoi)^deltai*exp(-deltai) > maxPeak) {
    mi <- maxPeak/c*(deltai/rhoi)^(-deltai)*exp(deltai)
  }
  initial.valuel <- c(sum(y), lambdai, deltai, rhoi) # m, lambda, delta, rho
  set <- which(is.na(fixed))
  initial.valuel <- initial.valuel[set]
  # Estimate parameters with MLE using the optim function
  modelTG1 <- optim(par = initial.valuel, fn = TiGo.lik, skew = 'L', maxPeak = maxPeak,error=error,
                    y = y, ultimateAdopt = ultimateAdopt,m = M, lambda=lambda, delta = delta, rho = rho,
                    method = method, control = control)
  initial <- modelTG1$par
  modelTG1 <- optim(par = initial, fn = TiGo.lik, skew = 'L', maxPeak = maxPeak,error=error,
                    y = y, ultimateAdopt = ultimateAdopt,m = M, lambda=lambda, delta = delta, rho = rho,
                    method = method, control = control)
  
  # Compare the likelihood of the right and left skewed models
  if (modelTG1$value < modelTG$value){
    modelTG <- modelTG1
  }
  if (!is.na(fixed[2]) && lambda < 0) modelTG <- modelTG1
  parameter <- fixed  
  parameter[set] <- modelTG$par  
  m <- parameter[1]
  lambda <- parameter[2]
  delta <- parameter[3]
  rho <- parameter[4]
  n <- length(y)
  T <- c(1 ): (n)
  I <- ifelse(lambda < 0, 1, 0)
  lowerGamma <- pgamma(rho, delta, lower = TRUE) * gamma(delta)
  gammaF <- gamma(delta) 
  c <- abs(lambda)*exp(delta*log(rho) - log(abs(lowerGamma - I * gammaF)))
  yf <- m * abs(lambda) * exp(-lambda*T*delta + delta*log(rho)-rho*exp(-lambda*T)-log(abs(lowerGamma - I * gammaF)))
  yf <- as.numeric(yf)
  # Calculate peak time and peak value
  peak <- m * exp(log(abs(lambda))- log(abs(lowerGamma - I * gammaF)) + delta * log(delta) - delta)
  peakTime <- -1/lambda*log(delta/rho)
  if (peakTime < 0) peakTime <- 0
  sigma <- sqrt(sum((log(yf)-log(y))^2)/length(y))
  out <- list(par = parameter, sigma = sigma, fit.y = as.numeric(yf), y = y, peakValue = as.numeric(peak), peakTime = peakTime, error=error)
  return(structure(out, class = "TiGo.model"))
}

print.TiGo.model <- function(model) {
  cat(paste("    lambda  = ", round(model$par[2], 4), "\n"))
  cat(paste("    delta  = ", round(model$par[3], 4), "\n"))
  cat(paste("    rho  = ", round(model$par[4], 4), "\n"))
  cat(paste("    m  = ", round(model$par[1], 4), "\n"))
  cat(paste("    peak time  = ", round(model$peakTime, 4), "\n"))
  cat(paste("    peak value  = ", round(model$peakValue, 4), "\n"))
}

pointForecast.TiGo <- function(model, h) {
  # -------------------------------------------------------------------
  # This function generates point forecasts from the TiGo model.
  # Args:
  #   model ---- the TiGo Diffusion model
  #   h ---- forecasting horizon
  # --------------------------------------------------------------------
  parameter <- model$par  
  error <- model$error
  m <- parameter[1]
  lambda <- parameter[2]
  delta <- parameter[3]
  rho <- parameter[4]
  n <- length(model$y)
  T <- (n+1) : (n + h)
  I <- ifelse(lambda < 0, 1, 0)
  lowerGamma <- pgamma(rho, delta, lower = TRUE) * gamma(delta)
  gammaF <- gamma(delta) 
  yf <- m * abs(lambda) * exp(-lambda*T*delta + delta*log(rho)-rho*exp(-lambda*T)-log(abs(lowerGamma - I * gammaF)))
  yf <- as.numeric(yf)
  y <- model$y
  out <- list(mean = yf, y = y, fitted = model$fit.y)
  return(structure(out, class = 'pointForecast.TiGo'))
}

forecast.TiGo <- function(model,h,level=c(80,95)) {
  # -------------------------------------------------------------------
  # This function generates forecasts from the TiGo model.
  # Args:
  #   model ---- the TiGo Diffusion model
  #   h ---- forecasting horizon
  #   level --- level of prediction intervals;default is 80% and 95%.
  #   error ---- error type; can be set to "normal" or "lognormal"
  # --------------------------------------------------------------------
  yf <- log(pointForecast.TiGo(model,h)$mean)
  error <- model$error
  sigma <- model$sigma
  interval <- matrix(0, h, 2*length(level))
  name <- matrix(0, 1, 2*length(level))
  for ( i in 1:length(level) ) {
    name[2*i] <- paste("Hi", level[i])
    name[2*(i-1)+1] <- paste("Lo", level[i])
  }
  # calculate quantiles according to the log-normal distribution.
  for (i in 1:length(level)){
    if (error == "normal") {
      interval[, 2*(i-1)+1] <-  qnorm(0.5 - level[i]/200, mean = yf, sd = sigma)
      interval[, 2*i] <- qnorm(0.5 + level[i]/200, mean = yf, sd = sigma)
    } else {
      interval[, 2*(i-1)+1] <- qlnorm(0.5 - level[i]/200, mean = yf, sd = sigma)
      interval[, 2*i] <- qlnorm(0.5 + level[i]/200, mean = yf, sd = sigma)
    }
  }
  colnames(interval) <- name
  out <- cbind(exp(yf), interval)
  result <- list(y = model$y, 'predictions' = out, median = exp(yf), level=level)
  return(structure(result, class = "TiGo.model")) 
}

plot.forecast.TiGo <- function(prediction,PI = TRUE,ylim = NULL,main=NULL,xlab=NULL,ylab=NULL){
  # ------------------------------------------------------------------------
  # This function creates a plot of the forecasts generated from a TiGo model.
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
      title <- 'Forecasts from TiGo'
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
      shadecols <- rev(colorspace::sequential_hcl(52 )[level - 49])
      if (h == 1){
        polygon(c(xxx + c(-0.5,0.5,0.5,-0.5)), c(rep(out1[c(j[(i-1)*2+1])],2),rep(out1[c(j[i*2])], 2)), 
                col  =  shadecols[i], border  =  FALSE) 
      } else {
        polygon(c(xxx, rev(xxx)), c(out1[, c(j[(i-1)*2+1])],rev(out1[, c(j[i*2])])), 
                col  =  shadecols[i], border  =  FALSE)
        lines(xxx,out[,1],col = 'blue',lwd = 2,xlab = xlab)
      }
    }
  }
  axis(1, pos=0) 
  axis(2, pos=0) 
}

