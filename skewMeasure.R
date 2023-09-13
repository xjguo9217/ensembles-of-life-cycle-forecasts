# -------------------------------
# Skewness measure based on mode
# -------------------------------

# This file contains R functions for calculating skewness values proposed in Arnold and Groeneveld (1995). They define the skewness of any distribution with a single-peaked density as Skew (F ) = 1 − 2F (t∗), where F (t∗) is the cdf evaluated at its mode t∗.

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

skewMeasure.lc <- function(model){
  if (class(model) == 'TiGo.model') {
    par <- model$par
    lambda <- par[2]
    delta <- par[3]
    rho <- par[4]
    if ((lambda > 0 & delta < rho) | (lambda < 0 & delta > rho)) {
      mode <- -1/lambda*log(delta/rho)
    } else {
      mode <- 0
    }
    if (!is.na(mode)&mode == 0) {
      skewMeasure = 1
    } else {
      if (is.na(mode)) {
        skewMeasure = -1
      } else {
      if (lambda > 0){
        F_mode <- (lower.inc.gamma(delta, rho)[1]-lower.inc.gamma(delta,rho*exp(-lambda*mode))[1])/(lower.inc.gamma(delta,rho)[1])
      }
      if (lambda < 0) {
        F_mode <- (lower.inc.gamma(delta,rho)[1]-lower.inc.gamma(delta,rho*exp(-lambda*mode))[1])/(lower.inc.gamma(delta,rho)[1]-gammaMpfr(delta))
      }
      skewMeasure <- as.numeric(1-2*F_mode)
      }
    }
    output <- list(mode = mode, skewness = skewMeasure, lcmodel = "TiGo")
    return(structure(output, class = 'skewMeasure'))
  }
  
  if (class(model) == "Bass.model"){
    M <- model$par[3]
    p <- model$par[1]
    q <- model$par[2]
    if (p < q){
      F_mode <- (1 - p/q)/2
      skewness = 1 - 2*F_mode
      output <- list(mode = mode, skewness = skewness, lcmodel = "Bass Diffusion Model", M = M)
      return(structure(output, class = 'skewMeasure'))
    }else{
      skewness =1
      output <- list(mode = 0, skewness = 1, lcmodel = "Bass Diffusion Model", M = M)
      return(structure(output, class = 'skewMeasure'))
    }
  }
  
  if (class(model) == "GSG"){
    parameter <- model$par
    M <- parameter[1]
    b <- parameter[2]
    beta <- 1/parameter[3]
    alpha <- parameter[4]
    A <- -beta*(alpha-1)^2
    B <- beta*alpha^2 + 3*alpha - 2
    C <- - 1/ beta - alpha
    D <- alpha*(-4 + 5*alpha - 4*beta + 4*alpha*beta + 2*alpha^2*beta + alpha^3*beta^2)
    if(is.na(D) | is.na(B) | D < 0 | B < 0){
      output <- list(mode = 0, skewness = NA, lcmodel = "GSG")
      #   if(D < 0 | B < 0 | (D > 0 & B > 0 & 2*A + B >0 & A + B + C < 0)){
      return(structure(output, class = 'skewMeasure'))
    }else{
      if (is.na(2*A + B) || is.na(A + B + C ) ||(2*A + B < 0 & A + B + C < 0)) {
        skewMeasure <- NA
      } else {
        xStar <- 2*C/(-B-sqrt(B^2 - 4*A*C))
        mode <- (-1/b)*log(xStar)
        mode2 <- mode*2
        F_mode <- (1-exp(-b*mode))/((1 + beta*exp(-b*mode))^alpha)
        skewMeasure <- 1- 2*F_mode
        if (is.na(mode) | mode <= 0 | mode == 0) skewMeasure <- NA
      }
      output <- list(mode = mode, skewness = skewMeasure, lcmodel = "GSG", M = M)
      return(structure(output, class = 'skewMeasure'))
    }
  }
  if (class(model) == "Trapezoid.model"){
    par <- model$par  
    a <- par[1]
    b <- par[2]
    c <- par[3]
    tau_1 <- par[4]
    tau_2 <- tau_1 + par[5]
    TT <- tau_2 - (a * tau_1 + b)/c
    
    Trap_cdf <- function(x) {
      if(0 <= x && x < tau_1) return(a * x^2 / 2 + b * x)
      if(tau_1 <= x && x < tau_2) return(a * tau_1^2 / 2 + b * tau_1 + (a * tau_1 + b) * (x - tau_1))
      if(tau_2 <= x && x < TT) return(- a * tau_1^2 / 2 + (a * tau_1 + b) * tau_2 + c * (x^2 - tau_2^2)/2 + (a * tau_1 + b - c * tau_2) * (x - tau_2))
      if(TT <= x)              return(- a * tau_1^2 / 2 + (a * tau_1 + b) * tau_2 + c * (TT^2 - tau_2^2)/2 + (a * tau_1 + b - c * tau_2)  * (TT - tau_2))
    }
    tmax <- tau_2 - (a*tau_1 + b)/c
    m <- Trap_cdf(TT)
    t_star <- (tau_1 + tau_2)/2
    skewMeasure <- 1 - 2 * Trap_cdf(t_star)/m
    output <- list(skewness = skewMeasure, lcmodel = "Trapezoid")
    return(structure(output, class = 'skewMeasure'))
  }
}

print.skewMeasure <- function(skewMeasure){
  if (skewMeasure$mode != 9999){
    cat(paste('Life Cycle model  :',skewMeasure$lcmodel,"\n"))
    cat(paste('Mode  :', round(skewMeasure$mode, 2),"\n"))
    cat(paste('Skewness  :', round(skewMeasure$skewness, 2),"\n"))
  }else{
    cat(paste('Life Cycle model  :',skewMeasure$lcmodel,"\n"))
    cat(paste('warning: the mode is at 0'))
  }
}