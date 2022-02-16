# ------------------------------------
# Tilted-Gompertz ETS model (TiGo-ETS)
# ------------------------------------

# This file contains R functions for fitting a TiGo-ETS model to data and generating forecasts from the model.  

# The gamma function in R base will return Inf when x in gamma(x) is greater than 171. Therefore, we need the Rmpfr function to calculate the gamma function and lower incomplete gamma function when x > 171.
library(Rmpfr)
gammaMpfr <- function(x){
  if (x > 170){
    return(gamma(as(x,"mpfr")))
  }else{
    return(gamma(x))
  }
}

lower.inc.gamma <- function(a, x) {
  lowerIncGamma <- pgamma(x, a, lower = TRUE) * gammaMpfr(a)[[1]]
  return(lowerIncGamma)
}

TiGo.ETS <- function(y = NULL, 
                     alpha = NULL,beta = NULL, 
                     initial.value = NULL, 
                     method = "Nelder-Mead") {
  # ------------------------------------------------------------------------------------------
  # This function fits the Tilted Gompertz Exponential Smoothing model.
  # Args:
  #   y ---- a numeric vector or time series     
  #   alpha ---- Value of alpha. If NULL, it is estimated.
  #   beta ---- Value of beta. If NULL, it is estimated.
  #   initial.value ---- initial value passed to optim
  #   method ---- method argument passed to optim
  #    An modelect of class "ets.TiGo".
  # ---------------------------------------------------------------------------------------------
  y1 <- y
  y <- log(y)
  
  # Fit a right-skewed TiGo model
  # Set initial values for MLE
  if(is.null(initial.value)){
    initial1 <-  TiGo.coefInt(y1, skew = 'R', alpha=alpha, beta=beta, initial.method = 'TiGo')
    initial <- matrix(initial1$coefs, 1, length(initial1$coefs))
  }else{
    initial <- initial.value
  }
  ets_lc <- optim(initial, TiGo.likelihoodFunc, y = y, 
                  skew = 'R', alpha = alpha, beta = beta,
                  method = method, control = list(maxit =2000))
  initial1 <- ets_lc$par
  ets_lc <- optim(initial1, TiGo.likelihoodFunc, y = y, 
                  skew = 'R', alpha = alpha, beta = beta,
                  method = method, control = list(maxit =2000))
  # If alpha or beta != 0, try another set of initial values for optim
  if ( is.null(alpha)||((!is.null(alpha) & alpha!=0)))  {
    initial <-  TiGo.coefInt(y1,skew = 'R', alpha=alpha, beta=beta)
    initial <- matrix(initial$coefs, 1, length(initial$coefs))
    ets_lc2 <- optim(initial, TiGo.likelihoodFunc, y = y, 
                     skew = 'R', alpha = alpha, beta = beta,
                     method = method, control = list(maxit =2000))
  }
  
  # Fit a left-skewed TiGo model
  if(is.null(initial.value)){
    initial1 <-  TiGo.coefInt(y1, skew = 'L', alpha=alpha, beta=beta, initial.method = 'TiGo')
    initial <- matrix(initial1$coefs, 1, length(initial1$coefs))
  }else{
    initial <- initial.value
  }
  ets_lc1 <- optim(initial, TiGo.likelihoodFunc, y = y, 
                   skew = 'L', alpha = alpha, beta = beta,
                   method = method, control = list(maxit =2000))
  initial1 <- ets_lc1$par
  ets_lc1 <- optim(initial1, TiGo.likelihoodFunc, y = y, 
                   skew = 'L', alpha = alpha, beta = beta,
                   method = method, control = list(maxit =2000))
  if ( is.null(alpha)||((!is.null(alpha) & alpha!=0)))  {
    initial <-  TiGo.coefInt(y1,skew = 'L',alpha=alpha, beta=beta)
    initial <- matrix(initial$coefs, 1, length(initial$coefs))
    ets_lc3 <- optim(initial, TiGo.likelihoodFunc, y = y, 
                     skew = 'L', alpha = alpha, beta = beta,
                     method = method, control = list(maxit =2000))
  }
  # Find the model with the smallest -loglikelihood.
  if ( is.null(alpha)||((!is.null(alpha) & alpha!=0)))   {
    likeli_values <- c(ets_lc$value,ets_lc1$value,ets_lc2$value,ets_lc3$value)
    nnn <- which.min(likeli_values)
    if (nnn == 2) ets_lc = ets_lc1
    if (nnn == 3) ets_lc = ets_lc2
    if (nnn == 4) ets_lc = ets_lc3
  } else {
    if (ets_lc1$value < ets_lc$value) ets_lc = ets_lc1
  }
  
  # Calculate in-sample one-step head forecast and the sd of errors
  para <- ets_lc$par
  if (!is.null(beta)) para <- cbind(beta, para)
  if (!is.null(alpha)) para <- cbind(alpha,para)
  fit <- Fit.TiGo.ETS(para, length(y), y) 
  fit.y <- fit$fit.y
  states <- cbind(fit$level, fit$trend)
  sigma <- sqrt(mean((y - fit.y)^2))
  para <- matrix(para, 1, length(para))
  colnames(para) <- c("alpha", "beta", "phi", "l0", "b0", "tau")
  
  # Calculate TiGo parameters lambda, delta, rho and m.
  phi <- para[3]
  tau <- para[6]
  lambda <- -log(phi)
  tt <- length(y)
  delta = log(tau)*phi/(log(phi)*(1-phi))
  states1 <- rbind(para[4:5],states)
  rho1 <- phi/(1-phi)*(states1[,2]-log(tau)/(1-phi)) # rho is time varying
  rho2 <- rho1[tt+1] # rho at the last time period
  lt <- states1[tt+1,1]
  bt <- states1[tt+1,2]
  logb0 = 1/(phi^tt)*(bt - (1-phi^(tt))/(1-phi)*log(tau))
  b0 = exp(logb0)
  l0 = lt - (phi - phi^(tt+1))/(1-phi)*log(b0) - log(tau)*phi/(1-phi)*(tt-1 - phi/(1 - phi) + phi^(tt)/(1-phi))
  l0 = exp(l0)
  TiGo_parameters = rbind(lambda, delta, matrix(rho1,length(rho1),1))
  c_0 <- NA
  lowerGamma <- lower.inc.gamma(delta,rho2)
  I = 1
  if (lambda > 0) {
    I = 0
  }
  gammaF <- gammaMpfr(delta)
  m_t <- log(states1[,1]) + rho2 + log((lowerGamma-I*gammaF)/lambda) - delta*log(rho2)
  m_t <- exp(m_t)
  m_t[is.na(m_t)] <- sum(exp(fit.y))
  
  # Predict peak time and peak value
  if (is.na(rho2)||rho2 < 0) {
    log_peak = 0
  } else {
    log_peak <- lt + rho2 - delta+ delta*(log(delta) - log(rho2))
  }
  peakValue <- exp(log_peak)
  rownamesRho <- NULL
  for (i in 1:dim(states1)[1]){
    rownamesRho <- rbind(rownamesRho, paste("rho", i-1))
  }
  rownames(TiGo_parameters) <- c("lambda", "delta", rownamesRho)
  if (is.na(rho2)|rho2 <0){
    timeToPeak <- 9999
  } else {
    timeToPeak <- -1/lambda*(log(delta) - log(rho2)) + length(y)
  }
  model_output <- list(y = y1, par = para, sigma = sigma, m_t = m_t, fit.y = exp(fit.y),states = states,
                       peakValue = peakValue, timeToPeak = timeToPeak, TiGo.parameters = TiGo_parameters)
  return(structure(model_output, class = "ets.TiGo"))
}

Fit.TiGo.ETS <- function(par, n, y) {
  # ------------------------------------------------------
  # This function calculates the in-sample fitted values.
  # -------------------------------------------------------
  alpha <- par[1];beta <- par[2];phi <- par[3]
  l0 <- par[4];b0 <- par[5];tau <- par[6]
  L <- matrix(0, n, 1);B <- matrix(0, n, 1)
  e <- matrix(0, n, 1);yf <- matrix(0, n, 1)
  yf[1] <- l0 + b0 * phi
  L[1] <- yf[1] + alpha * (y[1] - yf[1])
  B[1] <- b0 * phi + log(tau) + (beta) * (y[1] - yf[1])
  if (n > 1){
    for (i in 2:c(n)) {
      yf[i] <- L[i-1] + B[i-1] * phi
      L[i] <- yf[i] + alpha * (y[i] - yf[i])
      B[i] <- B[i-1] * phi + log(tau) + (beta)*(y[i] - yf[i])
    }
  }
  return(list(fit.y = yf, level = L, trend = B))
}

pointForecast.TiGo.ETS <- function(par, n) {
  # ------------------------------------------------------------------------------------------
  # This function calculates the fitted values given certain parameters.
  # ------------------------------------------------------------------------------------------
  phi <- par[3];l0 <- par[4];b0 <- par[5];tau <- par[6]
  L <- matrix(0, n, 1);B <- matrix(0, n, 1)
  e <- matrix(0, n, 1);yf <- matrix(0, n, 1)
  yf[1] <- l0 + b0*phi
  L[1] <- yf[1] 
  B[1] <- b0 * phi + log(tau)
  if (n > 1){
    for (i in 2:c(n)) {
      yf[i] <- L[i-1] + B[i-1] * phi
      L[i] <- yf[i]
      B[i] <- B[i-1] * phi + log(tau)
    }
  }
  return(list(fit.y = exp(yf), level = exp(L), trend = exp(B)))
}

forecast.TiGo.ETS <- function(model, h, PI = TRUE, level=c(80, 95)){
  # -------------------------------------------------------------------
  # This function generates forecasts from the Bass model.
  # Args:
  #   model ---- the Bass Diffusion model
  #   h ---- forecasting horizon
  #   PI ---- if FALSE only point forecasts will be generated
  #   level --- level of prediction intervals;default is 80% and 95%.
  # --------------------------------------------------------------------
  parameter <- model$par
  y <- model$y
  sigma <- model$sigma
  alpha <- model$par[1]
  beta <- model$par[2]
  phi <- c(model$par[3])
  tau <- model$par[6]
  if (is.null(y)) {
    l0 <- model$par[4]
    b0 <- model$par[5]
  } else {
    l0 <- model$states[nrow(model$states),1]
    b0 <- model$states[nrow(model$states),2]
    model$par[4:5] <-  model$states[nrow(model$states),]
  }
 
  if (PI == TRUE){
    interval <- matrix(0, h, 2*length(level))
    name <- matrix(0, 1, 2*length(level))
    for ( i in 1:length(level) ) {
      name[2*i] <- paste("Hi", level[i])
      name[2*(i-1)+1] <- paste("Lo", level[i])
    }
    var_t <- rep(NA, h)
    var_t1 <- rep(0, h)
    for (i in 1: (h-1)) {
      var_t1[i+1] <- var_t1[i] + (alpha + beta*(phi - phi^(i+1))/(1 - phi))^2
    }
    var_t <- (var_t1 + 1)*sigma^2
    forec1 <- pointForecast.TiGo.ETS(c(0,0,model$par[3:6]),h)
    median <- forec1$fit.y
    forec <- cbind(forec1$fit.y,forec1$level,forec1$trend)
    for (i in 1: length(level)) {
      p1 <- 0.5 - level[i]/200
      p2 <- 0.5 + level[i]/200
      interval[, 2*(i-1)+1] <- qlnorm(p1, log(median), sqrt(var_t))
      interval[, 2*i] <- qlnorm(p2,  log(median), sqrt(var_t)) 
    }
    colnames(interval) <- name
  }
  colnames(forec) <- c("point forecast", 'level', 'trend')
  
  if (PI == TRUE){
    out <- cbind(forec, interval)
  } else {
    out <- forec
  }
  result <- list(y = model$y, fit.y = model$fit.y, 'predictions' = out, median = median, mean = mean, level=level, PI = PI,var = var_t)
  return(structure(result, class = "forecast.Mlc")) 
}

TiGo.coefInt <- function(y, skew = 'R', alpha = NULL, beta = NULL, initial.method = 'lm'){
  # ------------------------------------------------------
  # This function initializes parameters passed to optim. 
  # -------------------------------------------------------
  fixed <- c(NA,NA)
  if (!is.null(alpha)) fixed[1] <- 1
  if (!is.null(beta)) fixed[2] <- 1
  checkbt <- 0
  n <- length(y)
  alpha <- ifelse(is.null(alpha),0.2,alpha)
  beta <- ifelse(is.null(beta),0.1,beta)
  phi <- ifelse(skew=='R',0.97,1.03)
  tau <- 0.95
  y1 <- log(y)
  # Use linear regression to find initial values for b0. 
  # The forecast function in R applies the same method.
  maxn <- min(10,n)
  if (maxn != 1) {
    fit <- lsfit(1:maxn,y1[1:maxn])
    l0 <- exp(fit$coef[1])
    b0 <- exp(fit$coef[2])
    if(abs(l0+b0) < 1e-8){
      l0 <- l0*(1+1e-3)
      b0 <- b0*(1-1e-3)
    }
  }
  # Estimate l0 using M. M is approximated by sum(y).
  if (initial.method == 'TiGo') {
    delta = log(tau)*phi/(log(phi)*(1-phi))
    rho = phi/(1-phi)*(log(b0) - log(tau)/(1-phi))
    lowerGamma <- lower.inc.gamma(delta,rho)#
    lambda <- -log(phi)
    I <- ifelse(lambda>0,0,1)
    gammaF <- gammaMpfr(delta)
    l0 <- exp(log(sum(y)) -( rho + log((lowerGamma-I*gammaF)/lambda) - delta*log(rho)))
  }
  lint <- l0
  bint <- b0
  
  # Check if rho < 0. Reset initial parameters if rho < 0 
  if (alpha == 0 && beta == 0) {
    bb <- bint
  } else {
    bb <- exp(Fit.TiGo.ETS(c(alpha, beta, phi, log(lint), log(bint), tau),length(y), y1)$trend[length(y)])
  }
  par <- check_b(c(alpha, beta, phi, bb, tau), skew)$par
  if (check_b(c(alpha, beta, phi, bb, tau), skew)$checkbt == 1) bint <- par[3]
  checkbt <- check_b(c(alpha, beta, phi, bb, tau), skew)$checkbt
  alpha <- par[1]
  beta <- par[2]
  
  coefs <- c(alpha, beta, phi, log(lint), log(bint), tau)
  if (length(which(fixed==1) > 0)) coefs <- coefs[-which(fixed==1)]
  output <- list(coefs = coefs, checkbt = checkbt)
  return(output)
}


TiGo.likelihoodFunc <- function(coef, y, skew='R', 
                                alpha = NULL,beta = NULL) {
  # -------------------------------------------------------------------
  # This function calculates the likelihood function. 
  # Args:
  #   coef --- parameters of a lc function (alpha, beta, phi, l0, b0, tau).
  #   y ---- time series
  #   skew --- skewness of the life cycle (R or L)
  #   alpha, beta --- fixed alpha and beta
  # --------------------------------------------------------------------
  alpha1 <- alpha; beta1 <- beta
  n <- length(y)
  #find if there are fixed parameters
  if (!is.null(beta)) coef <- c(beta,coef)
  if (!is.null(alpha)) coef <- c(alpha,coef)
  alpha <- coef[1]; beta  <- coef[2]; phi <- coef[3]
  l0    <- coef[4]; b0 <- coef[5]; tau   <- coef[6]
  
  # parameter constraints
  if (phi <= 0.0001) return(Inf)
  if (skew == 'R' & phi >= 0.9999)  return(Inf)
  if (skew == 'L' & phi <= 1.0001 ) return(Inf)
  if (tau >= 0.9999 ) return(Inf) 
  if (tau <= 0.0001 ) return(Inf)  
  if (beta >= 0.9999 | beta < 0 ) return(Inf)
  if (alpha >= 0.9999 | alpha < 0 ) return(Inf)
  
  # Calculate in-sample fitted value and likelihood
  fit <- Fit.TiGo.ETS(coef, n, y)
  yf <- fit$fit.y
  bt <- fit$trend[n]
  lt <- fit$level[n]
  lik <- log(sum((y - yf)^2)) * n
  bt <- exp(bt); b0 <- exp(b0); lt <- exp(lt); l0 <- exp(l0)
  if (is.na(lik) || is.na(bt) || bt == Inf || bt <= 0 || is.na(lt) || lt == Inf)
    return(Inf)
  
  if (is.null(alpha1)|| is.null(beta1)) {
    checkbt <- check_b(c(alpha, beta, phi, bt, tau), skew)$checkbt
    if (checkbt == 1) return(Inf)
  }
  if ((!is.null(alpha1)&&alpha1 ==0)&&(!is.null(beta1)&&beta1==0)) {
    checkbt <- check_b(c(alpha, beta, phi, b0, tau), skew)$checkbt
    if (checkbt == 1) return(Inf)
  }

  return(lik)
}

check_b <- function(par, skew) {
  # ------------------------------------------------------------------------------------------
  # This function checks if rho < 0. If rho < 0, reset initial values for alpha, beta and b0
  # Args:
  #   par ---- a numeric vector of parameters
  #   skew ---- skewness of the life cycle (R or L).
  # ------------------------------------------------------------------------------------------
  alpha <- par[1];beta <- par[2];phi <- par[3]
  b0 <- par[4];tau <- par[5]
  checkbt <- 0
  if (is.na(b0) || (phi > 1 & log(b0) - (log(tau)/(1-phi)) >= 0 )) {
    b0 <- tau^(1/(1-phi)) - 1e-5 
    alpha <- 0;beta <- 0
    checkbt <- 1
  }
  if (phi < 1 & log(b0) - (log(tau)/(1-phi)) <= 0){
    b0 <- tau^(1/(1-phi)) + 1e-5 
    alpha <- 0;beta <- 0
    checkbt <- 1
  }
  par <- c(alpha, beta, b0)
  output <- list(par = par, checkbt = checkbt)
  return (output)
} 

plot.forecast.Mlc <- function(prediction,PI = NULL,ylim = NULL,main=NULL,xlab=NULL,ylab=NULL){
  # --------------------------------------------------------------------------------
  # This function plots training data with point forecasts and prediction intervals. 
  #   prediction ---- prediction form the model
  #   PI ---- if FALSE only point forecasts will be plotted
  #   ylim, main, xlab, ylab --- graphic parameters
  # --------------------------------------------------------------------------------
  if (is.null(PI) == TRUE) PI <- prediction$PI
  level <- prediction$level
  y <- prediction$y
  h <- length(prediction$prediction[,1])
  out <- prediction$predictions
  if(any(is.na(out)))
    PI <- FALSE
  model_type <- prediction$model
  ytest <- prediction$ytest
  fit.y <- prediction$fit.y
  if (h == 1){
    out <- out[-c(2,3)]
  }else{
    if(PI == TRUE){
      out <- cbind(out[,1],out[,4:ncol(out)])
    } else {
      out <- matrix(out[,1])
    }
  }
  if (is.null(main)){
    title <- 'Forecasts from (M,Mlc,N)'
  }else{
    title <- main
  }
  if (is.null(ylim) ) {
    if (!is.null(ytest)){
      ylim <- c(min(min(out[!is.na(out)]), prediction$fit.y, y, ytest), 
                max(max(out[!is.na(out)]),prediction$fit.y, y, ytest))
    }else{
      ylim <- c(min(min(out[!is.na(out)]), y), max(max(out[!is.na(out)]), y))
    }
    if (ylim[2] == Inf){
      n0 <- which(out == Inf)
      ylim[2] <- max(out[-n0])
    }
  }
  plot(y,xlim = c(1,(length(y)+h)), ylim = ylim,type = 'l',main = title,xlab = xlab,ylab = ylab,axes = FALSE) 
  
  xxx <- (length(y)+1):(length(y)+h)
  if (h == 1){
    points(xxx,out[1],col = 'blue',lwd = 2)
  } else {
    lines(xxx,out[,1],col = 'black',lwd = 2,xlab = xlab, lty = 2)
  }
  if(PI == TRUE){
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
        points(xxx,out[1],col = 'blue',lwd = 2, pch=19)
      } else {
        polygon(c(xxx, rev(xxx)), c(out1[, c(j[(i-1)*2+1])],rev(out1[, c(j[i*2])])), 
                col  =  shadecols[i], border  =  FALSE)
        lines(xxx,out[,1],col = 'black',lwd = 2,xlab = xlab, lty = 2)
      }
    }
  }
  axis(1, pos=0) 
  axis(2, pos=0) 
}

