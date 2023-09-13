# ---------------------------------------------
# Forecasting Demand of Intel Microprocessors
# Using Four Time-invariant Diffusion Models
# (Multiple Model Types)
# ----------------------------------------------

# In this script, we generate Bayesian ensemble forecasts of Intel microprocessors sales.
# Candidate models considered here include the Bass, Gamma/shifted Gompertz, Trapezoid, and TiGo (time invariant) models.
library(zoo)
options(warn=-1)

# Load in the data of time series and start/end dates of each time series
data1 <- read.csv('Intel_actual_demand.csv')
startend_date <- read.csv('startend_date.csv')

timeSeriesIntel <- data1[,-c(1)] # remove the time columns
quant <- quant <- seq(0.01,0.99,0.01) # we consider 99 quantiles in model evaluation.

# Load in functions for the selected candidate model.
source('Source code/Bass.R')
source('Source code/GSG.R')
source('Source code/Trap.R')
source('Source code/TiGo-NLS.R')
alpha=0; beta=0

# Set up a function for calculating the pinball loss
pinball <- function(pred, act, q) {
  I = 1*(act>pred)
  l <- (act - pred)*q*I + (pred - act)*(1-q)*(1-I)
  return(l)
}

# Create a matrix of the number of time series for each time series at each time step ahead.
start <- startend_date$start
end <- startend_date$end

# Fit the selected model to all time series using the Maximum Likelihood Estimation method.
# For each of the 86 time series, then generate 187 steps ahead forecasts (the max length of the time series in the data set is 187) from time 0. Since these models are time invariant, 1-h step-ahead forecasts generated at time t are equal to t+1 - t+h step-ahead forecasts generated at time 0. 
candidateForecasts_bass <- matrix(NA,86, 187)
candidateForecasts_gsg <- matrix(NA,86, 187)
candidateForecasts_trap <- matrix(NA,86, 187)
candidateForecasts_tigo <- matrix(NA,86, 187)
sigma2Values_bass <- rep(NA,86)
sigma2Values_gsg <- rep(NA,86)
sigma2Values_trap <- rep(NA,86)
sigma2Values_tigo <- rep(NA,86)
min_value <-matrix(NA,86, 187)
max_value <- matrix(NA,86, 187)
n_drop <- c()
for (j in 1:86) {
  y1 <- as.numeric(timeSeriesIntel[,j])
  y <- c(na.omit(y1))
  y[which(y==0)] <- NA
  y <- (na.locf(y) + rev(na.locf(rev(y))))/2
  fitModel_bass <- Bass.fit(y)
  fitModel_gsg <- GSG.fit(y)
  fitModel_trap <- Trapezoid.fit(y)
  fitModel_tigo <- TiGo.NLS(y)
  sigma2Values_bass[j] <- (fitModel_bass$sigma)^2 # Store the variance of the errors to sigma2Values
  sigma2Values_gsg[j] <- (fitModel_gsg$sigma)^2
  sigma2Values_trap[j] <- (fitModel_trap$sigma)^2
  sigma2Values_tigo[j] <- (fitModel_tigo$sigma)^2
  # If the model suggests that the product has not passed its peak, we remove this product from the ensemble
  if (which.max(fitModel_bass$fit.y) == length(fitModel_bass$fit.y)) n_drop <- c(n_drop, j)
  # Generate forecasts
  fitModel_bass$y <- NULL
  fitModel_gsg$y <- NULL
  fitModel_trap$y <- NULL
  fitModel_tigo$y <- NULL
  forecasts_bass_j <- forecast.Bass(fitModel_bass,187,level = 99.9)$predictions
  forecasts_gsg_j <- forecast.GSG(fitModel_gsg,187,level = 99.9)$predictions
  forecasts_trap_j <- forecast.Trapezoid(fitModel_trap,187,level = 99.9)$predictions
  forecasts_tigo_j <- forecast.TiGo(fitModel_tigo,187,level = 99.9)$predictions
  candidateForecasts_bass[j,1:187] <- forecasts_bass_j[,1] # Store point forecasts to candidateForecasts
  candidateForecasts_gsg[j,1:187] <- forecasts_gsg_j[,1] 
  candidateForecasts_trap[j,1:187] <- forecasts_trap_j[,1] 
  candidateForecasts_tigo[j,1:187] <- forecasts_tigo_j[,1] 
  min_value[j,] <- min(forecasts_bass_j[,2],forecasts_gsg_j[,2],forecasts_trap_j[,2],forecasts_tigo_j[,2]) # Store the 0.1% and 99.9% quantiles 
  max_value[j,] <- max(forecasts_bass_j[,2],forecasts_gsg_j[,2],forecasts_trap_j[,2],forecasts_tigo_j[,2]) # Later we will use these values to calculate ensemble quantiles
}

# Set up accuracy matrix
accuracy_result <- array(NA,dim = c(86,187,17,99))
n_drop <- unique(n_drop)

# Calculate Bayesian ensembles and evaluate accuracy
for (i in 1:86) {
  y1 <- as.numeric(timeSeriesIntel[,i])
  y <- na.omit(y1)
  H1 <- length(y) - 8
  H2 <- length(y) - 17
  if (H1 < 0) next 
  # Generating rolling forecasts from t=0 to t=H1
  for (t in 0:H1) { 
    # Find products that are completed before t
    n_complete <- which(end < start[i]+t-1)
    n_complete <- n_complete[!n_complete%in%n_drop]
    if (length(n_complete) < 2) next # Move to the next one if there are less than two completed life cycles
    if (t > 0) {
      ytrain <- y[1:t]
    } else {
      ytain <- NULL
    }
    ytrain[which(ytrain==0)] <- NA
    ytrain <- (na.locf(ytrain) + rev(na.locf(rev(ytrain))))/2
    # Calculate weights on the completed life cycles
    weights <- rep(0,length(n_complete)*4)
    candidateForecasts1 <- rbind(candidateForecasts_bass[n_complete,],candidateForecasts_gsg[n_complete,],candidateForecasts_trap[n_complete,],candidateForecasts_tigo[n_complete,])
    sigma2Values1 <- c(sigma2Values_bass[n_complete],sigma2Values_gsg[n_complete],sigma2Values_trap[n_complete],sigma2Values_tigo[n_complete])
    # The following loop calculates the likelihood of the on-going lifecycle’s past data according to candidate model mm
    for (mm in 1:(length(n_complete)*4)) {
      if (t>0) {
        weights[mm] <- sum(log(dlnorm(ytrain,mean = log(candidateForecasts1[mm,1:t]), sd = sqrt(sigma2Values1[mm])))) 
      } else {
        weights[mm] <- 1
      }
    }
    if (sum(exp(weights)) == 0) {
      weights <- exp(weights-max(weights))/sum(exp(weights-max(weights)))
    } else {
      weights <- exp(weights)/sum(exp(weights))
    }
    
    # Generate forecasts for 1% to 99% quantiles
    if (t <= H2)  ntest <- 17
    if (t > H2 & t <= H1)  ntest <- 8
    forecast_matrix <- candidateForecasts1[,(t+1):(t+ntest)]
    fit_matrix <- as.matrix(candidateForecasts1[,1:t])
    forecast_ntest <- matrix(NA,ntest,99)
    # This is the function for calculating the ensemble CDF
    cdfFunc <- function(x,mean_value,sigma_value) {
      return_cdf <- 0
      for (k in 1:(length(n_complete)*4)) {
        return_cdf <- return_cdf + weights[k]*plnorm(x,log(mean_value[k]),sigma_value[k])
      }
      return(return_cdf)
    }
    # We set up a vector of x_values and calculate each element's ensemble CDF.
    # For each of the 1-99% quantiles p, find the nearest CDF calculated in the previous step. Its corresponding x value is the p quantile.
    min_value1 <- min_value[n_complete,(t+1):(t+ntest)]
    max_value1 <- max_value[n_complete,(t+1):(t+ntest)]
    for (kk in 1:ntest) {
      increment <- (max(max_value1[,kk]) - min(min_value1[,kk]))/2000
      x_values <- seq(min(min_value1[,kk]),max(max_value1[,kk]),increment)
      cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(sigma2Values1))
      # Check if the values in x_value can roughly cover the entire distribution. 
      # If not, reset x_values and repeat previous steps.
      if (length(unique(cdf_values)) <= 10) {
        x_values <- seq(0,1000,0.1)
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(sigma2Values1))
      }
      cdf_1_n <- which(cdf_values >= 0.5)[1]
      nn = 0
      while(!is.na(cdf_1_n) & cdf_1_n < 500) { # This is the case when half of the values fall above the 99% quantile.
        nn = nn+1
        increment <- (x_values[cdf_1_n] - min(min_value1[,kk]))/2000
        if (increment == 0) break
        x_values <- seq(min(min_value1[,kk]),x_values[cdf_1_n],increment)
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(sigma2Values1))
        cdf_1_n <- which(cdf_values >= 0.5)[1]
      }
      nn = 0
      while (round(max(cdf_values),2) < 0.99) { # This is the case when the max value is below the 99% quantile
        nn = nn+1
        increment <- x_values[which.max(cdf_values)] - x_values[which.min(abs(cdf_values - max(cdf_values)+0.01))]
        if (sum(weights*max_value1[,kk])+increment*2000 == Inf|nn > 10) break # stop after running the loop for 10 times
        x_values  <- c(x_values, seq(tail(x_values,1),tail(x_values,1)+increment*2000,increment))
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(sigma2Values1))
      }
      # Calculate quantiles
      for (p in 1:99)
        forecast_ntest[kk,p] <- x_values[which.min(abs(cdf_values - p/100))]
    }
    
    # Calculate pinball losses
    ytest <- y[(t+1):(t+ntest)]
    for (mm in 1:99) {
      accuracy_result[i,t+1,1:ntest,mm] <- pinball(forecast_ntest[,mm], ytest, mm/100)
    }
  }
}

# Print out average accuracy results
result_vec <- t(2*c(mean(rowMeans(apply(accuracy_result[,,1:8,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                    mean(rowMeans(apply(apply(accuracy_result[,,1:8,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                    mean(rowMeans(apply(accuracy_result[,,9:17,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                    mean(rowMeans(apply(apply(accuracy_result[,,9:17,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE)))
colnames(result_vec) <- c("1-8 steps MAE","1-8 steps MCRPS","9-17 steps MAE","13-24 steps MCRPS")
result_vec