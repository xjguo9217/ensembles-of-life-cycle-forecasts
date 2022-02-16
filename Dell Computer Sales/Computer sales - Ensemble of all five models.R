# -----------------------------------------------------------
# Forecasting Dell Computer Sales Using Bayesian Ensemble 
# of All Five Model Types (Including the TiGo-ETS model)
# -----------------------------------------------------------

# In this script, we generate Bayesian ensemble forecasts of search interests of social networking websites.
# Candidate models considered here include the Bass, Gamma/shifted Gompertz, Trapezoid, TiGo (time invariant) and TiGo-ETS models.

# Load in the lubridate package
library(lubridate)
options(warn=-1)

# Load in data and unnormalize the normalized life cycles
data1 <- read.csv('dell_data_truncated.csv')
weekly_order <- read.csv('weekly_M.csv')
timeSeriesDell <- data1[,-c(1:3)]
for (i in 1:170) {
  y1 <- as.numeric(timeSeriesDell[i, ])
  timeSeriesDell[i, ] <- as.numeric(timeSeriesDell[i, ])*weekly_order[i,1]
}

quant <- quant <- seq(0.01,0.99,0.01) # we consider 99 quantiles in model evaluation.

# Load in functions for the selected candidate model.
source('Bass.R')
source('GSG.R')
source('Trap.R')
source('TiGo-ETS.R')
alpha=0.3; beta=0.0003

# Set up a function for calculating the pinball loss
pinball <- function(pred, act, q) {
  I = 1*(act>pred)
  l <- (act - pred)*q*I + (pred - act)*(1-q)*(1-I)
  return(l)
}

# Convert products' start and end dates into objects of class "POSIXlt"
# Later we will use these dates to determine completed life cycles.
start <- strptime(data1[,1],"%Y-%m-%d",tz = 'GMT')
end <- strptime(data1[,2],"%Y-%m-%d",tz = 'GMT')

# Fit the selected model to all time series using Maximum Likelihood Estimation.
# For each of the 170 time series, then generate 82 steps ahead forecasts (the max length of the time series in the data set is 82) from the four time invariant models from time 0. For the TiGo-ETS model, we will update states and generate forecasts later.
candidateForecasts_bass <- matrix(NA,170, 82)
candidateForecasts_gsg <- matrix(NA,170, 82)
candidateForecasts_trap <- matrix(NA,170, 82)
candidateForecasts_tigo <- matrix(NA,170, 82)
candidateForecasts_ets<- matrix(NA,170, 82)
sigma2Values_bass <- rep(NA,170)
sigma2Values_gsg <- rep(NA,170)
sigma2Values_trap <- rep(NA,170)
sigma2Values_tigo <- rep(NA,170)
sigma2Values_ets <- rep(NA,170)
min_value <-matrix(NA,170, 82)
max_value <- matrix(NA,170, 82)
par_vec <- matrix(NA,170,4)
n_drop <- c()
for (j in 1:170) {
  y1 <- as.numeric(timeSeriesDell[j,])
  y <- c(na.omit(y1))
  y[y<=0] <- min(y[y>0])/2
  fitModel_bass <- Bass.fit(y)
  fitModel_gsg <- GSG.fit(y)
  fitModel_trap <- Trapezoid.fit(y)
  fitModel_tigo <- TiGo.ETS(y,alpha=0,beta=0)
  fitModel_ets <- TiGo.ETS(y,alpha=alpha,beta=beta)
  par_vec[j,] <- fitModel_ets$par[3:6] # store parameters
  sigma2Values_bass[j] <- (fitModel_bass$sigma)^2 # Store the variance of the errors to sigma2Values
  sigma2Values_gsg[j] <- (fitModel_gsg$sigma)^2
  sigma2Values_trap[j] <- (fitModel_trap$sigma)^2
  sigma2Values_tigo[j] <- (fitModel_tigo$sigma)^2
  sigma2Values_ets[j] <- (fitModel_ets$sigma)^2
  # If the model suggests that the product has not passed its peak, we remove this product from the ensemble
  if (tail(exp(fitModel_ets$states[,2]),1)>1) n_drop <- c(n_drop,j)
  # Generate forecasts
  fitModel_bass$y <- NULL
  fitModel_gsg$y <- NULL
  fitModel_trap$y <- NULL
  fitModel_tigo$y <- NULL
  forecasts_bass_j <- forecast.Bass(fitModel_bass,82,level = 99.9)$predictions
  forecasts_gsg_j <- forecast.GSG(fitModel_gsg,82,level = 99.9)$predictions
  forecasts_trap_j <- forecast.Trapezoid(fitModel_trap,82,level = 99.9)$predictions
  forecasts_tigo_j <- forecast.TiGo.ETS(fitModel_tigo,82,level = 99.9)$predictions[,-c(2,3)]
  candidateForecasts_bass[j,1:82] <- forecasts_bass_j[,1] # Store point forecasts to candidateForecasts
  candidateForecasts_gsg[j,1:82] <- forecasts_gsg_j[,1] 
  candidateForecasts_trap[j,1:82] <- forecasts_trap_j[,1] 
  candidateForecasts_tigo[j,1:82] <- forecasts_tigo_j[,1] 
  min_value[j,] <- min(forecasts_bass_j[,2],forecasts_gsg_j[,2],forecasts_trap_j[,2],forecasts_tigo_j[,2]) # Store the 0.1% and 99.9% quantiles 
  max_value[j,] <- max(forecasts_bass_j[,2],forecasts_gsg_j[,2],forecasts_trap_j[,2],forecasts_tigo_j[,2]) # Later we will use these values to calculate ensemble quantiles
}
TiGoModel1 <- fitModel_ets

# Set up the accuracy matrix
accuracy_result <- array(NA,dim = c(170,82,24,99))
n_drop <- unique(n_drop)

# Calculate Bayesian ensembles and evaluate accuracy
for (i in 1:170) {
  y1 <- as.numeric(timeSeriesDell[i,])
  y <- na.omit(y1)
  H1 <- length(y) - 8
  H2 <- length(y) - 17
  if (H1 < 0) next 
  # Generating rolling forecasts from t=0 to t=H1
  for (t in 0:H1) { 
    if (t <= H2)  ntest <- 17
    if (t > H2 & t <= H1)  ntest <- 8
    # Find products that are completed before t
    n_complete <- which((end < start[i]+days(7*(t-1))))
    n_complete <- n_complete[!n_complete%in%n_drop]
    if (length(n_complete) < 2) next # Move to the next one if there are less than two completed life cycles
    if (t > 0) {
      ytrain <- y[1:t]
      # We replace zeros in these series with the minimum positive value in the series, divided by 2.
      # If all the data points up to t are zeros, we replace them with the minimum positive value in the completed life cycles, divided by 2.
      if (any(ytrain<=0)) {
        if (all(ytrain == 0)) {
          prior_curves <- timeSeriesDell[n_complete,]
          prior_curves[prior_curves == 0] <- 1e10
          ytrain[which(ytrain==0)] <- min(prior_curves,na.rm=TRUE)/2
        } else {
          ytrain[which(ytrain==0)] <- min(ytrain[ytrain>0])/2
        }
      }
    }
    
    update_median <- matrix(NA,length(n_complete),t)
    sigma2Values_TiGo1 <- rep(NA, length(n_complete))
    min_value2 <-matrix(NA,length(n_complete), ntest)
    max_value2 <- matrix(NA,length(n_complete), ntest)
    sigma2Values_pred <- matrix(NA, length(n_complete),ntest)
    forecast_matrix <- matrix(NA, length(n_complete),ntest)
    n_remove <- c()
    # For each candidate model kk, update its states bt and lt using the available data from the current time series.
    for (kk in 1:length(n_complete)) {
      par <- c(alpha,beta,par_vec[n_complete[kk],])
      TiGoModel1$par <- par
      # update states
      if (t > 0) {
        fit <- Fit.TiGo.ETS(par, length(ytrain), log(ytrain))
        fitted <- fit$fit.y
        TiGoModel1$states <- cbind((fit$level),(fit$trend)) 
        TiGoModel1$sigma <- sigma2Values_tigo[n_complete[kk]]
        TiGoModel1$y <- ytrain
      } else {
        TiGoModel1$y <- NULL
      }
      if (t == 0 && (par[3] >= 1 & par[3]/(1-par[3])*(par[5]-log(par[6])/(1-par[3]))<=0)) n_remove <- c(n_remove,kk)
      if (t > 0 && (par[3] >= 1 & par[3]/(1-par[3])*(tail(fit$trend,1)-log(par[6])/(1-par[3]))<=0)) n_remove <- c(n_remove,kk) # remove the model if rho_t < 0 and therefore it is not a life cycle model
      # Find the median and variance used to calculate the likelihood of the on-going life cycle's past data according to candidate model kk.
      if (t>0) update_median[kk,1:t] <- exp(fitted) 
      sigma2Values_TiGo1[kk] <- sigma2Values_tigo[n_complete[kk]] 
      forc1 <- forecast.TiGo.ETS(TiGoModel1,h=ntest,level = 99.9)
      min_value2[kk,] <- forc1$predictions[,4] # Store the 0.1% and 99.9% quantiles
      max_value2[kk,] <- forc1$predictions[,5] # Later we will use these values to calculate ensemble quantiles
      forecast_matrix[kk,] <- forc1$predictions[,1] # Store point forecasts to forecast_matrix
      sigma2Values_pred[kk,] <- forc1$var # Store the variance for the 1-h steps ahead forecasts
    }
    # Calculate weights on the completed life cycles
    weights <- rep(0,length(n_complete)*5)
    if (t>0) candidateForecasts1 <- rbind(as.matrix(candidateForecasts_bass[n_complete,1:t]),as.matrix(candidateForecasts_gsg[n_complete,1:t]),as.matrix(candidateForecasts_trap[n_complete,1:t]),as.matrix(candidateForecasts_tigo[n_complete,1:t]),as.matrix(update_median))
    candidateForecasts2 <- rbind(candidateForecasts_bass[n_complete,],candidateForecasts_gsg[n_complete,],candidateForecasts_trap[n_complete,],candidateForecasts_tigo[n_complete,])
    sigma2Values1 <- c(sigma2Values_bass[n_complete],sigma2Values_gsg[n_complete],sigma2Values_trap[n_complete],sigma2Values_tigo[n_complete],sigma2Values_TiGo1)
    sigma2Values1_pred <- c(sigma2Values_bass[n_complete],sigma2Values_gsg[n_complete],sigma2Values_trap[n_complete],sigma2Values_tigo[n_complete])
    # The following loop calculates the likelihood of the on-going lifecycle’s past data according to candidate model mm
    for (mm in 1:(length(n_complete)*5)) {
      if (t>0) {
        weights[mm] <- sum(log(dlnorm(ytrain,mean = log(candidateForecasts1[mm,1:t]), sd = sqrt(sigma2Values1[mm])))) 
      } else {
        weights[mm] <- 1
      }
    }
    weights[((length(n_complete)*4 +1):(length(n_complete)*5))[n_remove]] <- -Inf
    if (sum(exp(weights)) == 0) {
      weights <- exp(weights-max(weights))/sum(exp(weights-max(weights)))
    } else {
      weights <- exp(weights)/sum(exp(weights))
    }
    
    # Generate forecasts for 1% to 99% quantiles
    forecast_matrix <- rbind(candidateForecasts2[,(t+1):(t+ntest)],forecast_matrix)
    forecast_ntest <- matrix(NA,ntest,99)
    # This is the function for calculating the ensemble CDF
    cdfFunc <- function(x,mean_value,sigma_value) {
      return_cdf <- 0
      for (k in 1:(length(n_complete)*5)) {
        return_cdf <- return_cdf + weights[k]*plnorm(x,log(mean_value[k]),sigma_value[k])
      }
      return(return_cdf)
    }
    # We set up a vector of x_values and calculate each element's ensemble CDF.
    # For each of the 1-99% quantiles p, find the nearest CDF calculated in the previous step. Its corresponding x value is the p quantile.
    max_value2[max_value2==Inf] <- 1e5
    min_value2[min_value2==Inf] <- 1e5
    min_value1 <- rbind(min_value[n_complete,(t+1):(t+ntest)],min_value2)
    max_value1 <- rbind(max_value[n_complete,(t+1):(t+ntest)],max_value2)
    for (kk in 1:ntest) {
      increment <- (sum(weights*max_value1[,kk]) - sum(weights*min_value1[,kk]))/2000
      x_values <- seq(sum(weights*min_value1[,kk]),sum(weights*max_value1[,kk]),increment)
      cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(c(sigma2Values1_pred,sigma2Values_pred[,kk])))
      # Check if the values in x_value can roughly cover the entire distribution. 
      # If not, reset x_values and repeat previous steps.
      if (length(unique(cdf_values)) <= 10) {
        x_values <- seq(0,1000,0.1)
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(c(sigma2Values1_pred,sigma2Values_pred[,kk])))
      }
      cdf_1_n <- which(cdf_values >= 0.5)[1]
      nn = 0
      while(!is.na(cdf_1_n) & cdf_1_n < 500) { # This is the case when half of the values fall above the 99% quantile.
        nn = nn+1
        increment <- (x_values[cdf_1_n] - min(min_value1[,kk]))/2000
        if (increment == 0) break
        x_values <- seq(min(min_value1[,kk]),x_values[cdf_1_n],increment)
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(c(sigma2Values1_pred,sigma2Values_pred[,kk])))
        cdf_1_n <- which(cdf_values >= 0.5)[1]
      }
      nn = 0
      while (round(max(cdf_values),2) < 0.99) { # This is the case when the max value is below the 99% quantile
        nn = nn+1
        increment <- x_values[which.max(cdf_values)] - x_values[which.min(abs(cdf_values - max(cdf_values)+0.01))]
        if (sum(weights*max_value1[,kk])+increment*2000 == Inf|nn > 10) break # stop after running the loop for 10 times
        x_values  <- c(x_values, seq(tail(x_values,1),tail(x_values,1)+increment*2000,increment))
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(c(sigma2Values1_pred,sigma2Values_pred[,kk])))
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
result_vec <- t(as.matrix(2*c(mean(rowMeans(apply(accuracy_result[,,1:8,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                              mean(rowMeans(apply(apply(accuracy_result[,,1:8,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                              mean(rowMeans(apply(accuracy_result[,,9:17,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                              mean(rowMeans(apply(apply(accuracy_result[,,9:17,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE))))
colnames(result_vec) <- c('1-8 steps MAE','1-8 steps CRPS','9-17 steps MAE','9-17 steps CRPS')
result_vec
