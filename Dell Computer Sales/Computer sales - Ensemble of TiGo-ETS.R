# --------------------------------------------------------------------------
# Forecasting Computer Sales Using the Bayesian Ensemble of TiGo-ETS Models
# (Single Model Type)
# --------------------------------------------------------------------------

# In this script, we generate Bayesian ensemble forecasts of Dell computer sales.
# We consider TiGo-ETS models as candidate models.

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

nDell <- nrow(timeSeriesDell)
quant <- quant <- seq(0.01,0.99,0.01) # we consider 99 quantiles in model evaluation.

# Load in functions for TiGo models.
source('TiGo-ETS.R')

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

# Fit the model to all time series using the Maximum Likelihood Estimation approach.
sigma2Values_TiGo <- rep(NA,170)
n_drop <- c()
par_vec <- matrix(NA,170,4)
# Set alpha and beta
alpha <- 0.3
beta <- 0.0003
for (j in 1:170) {
  y1 <- as.numeric(timeSeriesDell[j,])
  y <- c(na.omit(y1))
  y[y<=0] <- min(y[y>0])/2
  TiGoModel <- TiGo.ETS(y,alpha=alpha,beta=beta)
  par_vec[j,] <- TiGoModel$par[3:6] # store parameters
  sigma2Values_TiGo[j] <- (TiGoModel$sigma)^2 # Store the variance of the errors to sigma2Values
  # If the model suggests that the product has not passed its peak, we remove this product from the ensemble
  if (tail(exp(TiGoModel$states[,2]),1)>1) n_drop <- c(n_drop,j)
}

# Set up accuracy matrix
accuracy_result <- array(NA,dim = c(170,82,34,99))
n_drop <- unique(n_drop)
TiGoModel1 <- TiGoModel

# Calculate Bayesian ensembles and evaluate the accuracy
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
    ytrain <- y[1:t]
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
    
    # Calculate weights on the completed life cycles
    weights <- rep(0,length(n_complete))
    update_median <- matrix(NA,length(n_complete),t)
    sigma2Values_TiGo1 <- rep(NA, length(n_complete))
    min_value1 <-matrix(NA,length(n_complete), ntest)
    max_value1 <- matrix(NA,length(n_complete), ntest)
    sigma2Values_pred <- matrix(NA, length(n_complete),ntest)
    forecast_matrix <- matrix(NA, length(n_complete),ntest)
    n_remove <- c()
    # For each candidate model kk, update its states bt and lt using the available data from the current time series.
    for (kk in 1:length(n_complete)) {
      par <- c(alpha,beta,par_vec[n_complete[kk],])
      TiGoModel1$par <- par
      # update states
      fit <- Fit.TiGo.ETS(par, length(ytrain), log(ytrain))
      fitted <- fit$fit.y
      TiGoModel1$states <- cbind((fit$level),(fit$trend)) 
      TiGoModel1$sigma <- sigma2Values_TiGo[n_complete[kk]]
      TiGoModel1$y <- ytrain
      if (t == 0) TiGoModel1$y <- NULL
      if (par[3] >= 1 & par[3]/(1-par[3])*(tail(fit$trend,1)-log(par[6])/(1-par[3]))<=0) n_remove <- c(n_remove,kk) # remove the model if rho_t < 0 and therefore it is not a life cycle model
      # Find the median and variance vectors used to calculate the likelihood of the on-going life cycle's past data according to candidate model kk.
      if (t>0) update_median[kk,1:t] <- exp(fitted) 
      sigma2Values_TiGo1[kk] <- sigma2Values_TiGo[n_complete[kk]] 
      forc1 <- forecast.TiGo.ETS(TiGoModel1,h=ntest,level = 99.9)
      min_value1[kk,] <- forc1$predictions[,4] # Store the 0.1% and 99.9% quantiles
      max_value1[kk,] <- forc1$predictions[,5] # Later we will use these values to calculate ensemble quantiles
      forecast_matrix[kk,] <- forc1$predictions[,1] # Store point forecasts to forecast_matrix
      sigma2Values_pred[kk,] <- forc1$var # Store the variance for the 1-h steps ahead forecasts
    }
    
    # The following loop calculates the likelihood of the on-going lifecycle’s past data according to candidate model mm
    for (mm in 1:length(n_complete)) {
      if (t>0) {
        weights[mm] <- sum(log(dlnorm(ytrain,mean = log(update_median[mm,1:t]), sd = sqrt(sigma2Values_TiGo1[mm]))))
      } else {
        weights[mm] <- 1
      }
    }
    weights[n_remove] <- -Inf
    if (sum(exp(weights)) == 0) {
      weights <- exp(weights-max(weights))/sum(exp(weights-max(weights)))
    } else {
      weights <- exp(weights)/sum(exp(weights))
    } 
    
    # Generate forecasts for 1% to 99% quantiles
    ytest <- y[(t+1):(t+ntest)]
    forecast_ntest <- matrix(NA,ntest,99)
    
    # Replace Inf 99.9 quantiles with a large number.   
    max_value1[max_value1==Inf] <- 1e5
    min_value1[min_value1==Inf] <- 1e5
    # This is the function for calculating the ensemble CDF
    cdfFunc <- function(x,mean_value,sigma_value) {
      return_cdf <- 0
      for (k in 1:length(n_complete)) {
        return_cdf <- return_cdf + weights[k]*plnorm(x,log(mean_value[k]),sigma_value[k])
      }
      return(return_cdf)
    }
    # We set up a vector of x_values and calculate each element's ensemble CDF.
    # For each of the 1-99% quantiles p, find the nearest CDF calculated in the previous step. Its corresponding x value is the p quantile.
    for (kk in 1:ntest) {
      increment <- (sum(weights*max_value1[,kk]) - sum(weights*min_value1[,kk]))/2000
      x_values <- seq(sum(weights*min_value1[,kk]),sum(weights*max_value1[,kk]),increment)
      cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(sigma2Values_pred[,kk]))
      # Check if the values in x_value can roughly cover the entire distribution. 
      # If not, reset x_values and repeat previous steps.
      if (length(unique(cdf_values)) <= 10) {
        x_values <- seq(0,1000,0.1)
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(sigma2Values_pred[,kk]))
      }
      cdf_1_n <- which(cdf_values >= 0.5)[1]
      nn = 0
      while(!is.na(cdf_1_n) & cdf_1_n < 500) { # This is the case when half of the values fall above the 99% quantile.
        nn = nn+1
        increment <- (x_values[cdf_1_n] - min(min_value1[,kk]))/2000
        if (increment == 0) break
        x_values <- seq(min(min_value1[,kk]),x_values[cdf_1_n],increment)
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(sigma2Values_pred[,kk]))
        cdf_1_n <- which(cdf_values >= 0.5)[1]
      }
      nn = 0
      while (round(max(cdf_values),2) < 0.99) { # This is the case when the max value is below the 99% quantile
        nn = nn+1
        increment <- x_values[which.max(cdf_values)] - x_values[which.min(abs(cdf_values - max(cdf_values)+0.01))]
        if (sum(weights*max_value1[,kk])+increment*2000 == Inf|nn > 10) break # stop after running the loop for 10 times
        x_values  <- c(x_values, seq(tail(x_values,1),tail(x_values,1)+increment*2000,increment))
        cdf_values <- cdfFunc(x_values,forecast_matrix[,kk],sqrt(sigma2Values_pred[,kk]))
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
colnames(result_vec) <- c('1-8 steps MAE','1-8 steps CRPS','9-17 steps MAE','9-17 steps CRPS')
result_vec

# split the results for training and testing set
start_allts <- min(start)
end_allts <- max(end)
n_train <- round(interval(start_allts, end_allts) %/% weeks(1)/2)
end_train <- sort(start)[170/2]

end_train_ts <- interval(start, end_train) %/% weeks(1)+1

accuracy_train <- array(NA,dim = c(170,82,17,99))
accuracy_test <- array(NA,dim = c(170,82,17,99))
for (i in 1:170) {
  this_end_train <- end_train_ts[i]
  if (this_end_train <= 1) {
    accuracy_test[i,,1:17,] <- accuracy_result[i,,1:17,] 
  } else {
    if (this_end_train > 82) {
      accuracy_train[i,,1:17,] <- accuracy_train[i,,1:17,] 
    } else {
      accuracy_train[i,1:(this_end_train-1),1:17,] <- accuracy_result[i,1:(this_end_train-1),1:17,] 
      accuracy_test[i,this_end_train:82,1:17,] <- accuracy_result[i,this_end_train:82,1:17,] 
    }
  }
}

result_test <- t(2*c(mean(rowMeans(apply(accuracy_test[,,1:8,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                   mean(rowMeans(apply(apply(accuracy_test[,,1:8,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                   mean(rowMeans(apply(accuracy_test[,,9:17,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                   mean(rowMeans(apply(apply(accuracy_test[,,9:17,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE)))
colnames(result_test) <- c("1-8 steps MAE","1-8 steps MCRPS","9-17 steps MAE","9-17 steps MCRPS")
result_test
