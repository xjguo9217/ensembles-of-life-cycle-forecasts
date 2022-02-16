# --------------------------------------------------------------
# Forecasting Search Interests of Social Networking Websites 
# Using the TiGo Model Updated with Kalman Filter
# --------------------------------------------------------------
# In this script, we generate forecasts for search interest of social networking websites 
# using the TiGo-Kalman filter method.

# Load in the lubridate package
library(lubridate)

# Load in the data of time series and start/end dates of each time series
data1 <- read.csv('GoogleTrendsTimeSeriesSocial-May-2021.csv')
startend_dates <- read.csv('WikiListSocial-startenddays-May-2021.csv') 

timeSeriesSearchInterests <- data1[,-c(1:2)] # remove the time columns
quant <- quant <- seq(0.01,0.99,0.01) # we consider 99 quantiles in model evaluation.

# Load in functions for TiGo models.
source('TiGo-ETS.R')

# Set up a function for calculating the pinball loss
pinball <- function(pred, act, q) {
  I = 1*(act>pred)
  l <- (act - pred)*q*I + (pred - act)*(1-q)*(1-I)
  return(l)
}

# Create a matrix of the number of time series for each time series at each time step ahead.
start <- as.Date('2004-01-01')
for (i in 1:161) {
  y1 <- as.numeric(timeSeriesSearchInterests[,i])
  ll <- as.numeric(attr(na.omit(y1),"na.action"))
  if (length(ll) == 1 && ll == 1) {
    ll <- 1
  } else {
    if (any(diff(ll)!=1)) {
      ll <- which(diff(ll)!=1)[1]
    } else {
      ll <- tail(ll,1)
    }
  }
  if (length(ll) == 0) ll <- 0
  start <- c(start,ymd(as.Date('2004-01-01')) %m+% months(ll))
}

start <- start[-1]

end <-  as.Date(strptime(startend_dates$end_date,"%d/%m/%y",tz = 'GMT')) # for less than 5
end[which(is.na(end))] <- as.Date('2021-06-01')

# Pre-process data. Remove search interests data for periods that are after the 5 months + end date listed on Wikipedia. 
for (i in 1:161) {
  y1 <- timeSeriesSearchInterests[,i]
  if (end[i] != '2021-06-01') {
    n_end <- which(data1$month == startend_dates[i,7]&data1$year == startend_dates[i,8])
    if (n_end+4 <= 209) timeSeriesSearchInterests[(n_end+4):209,i] <- NA # for less than 5
  }
}

# Fit the model to all time series using the Maximum Likelihood Estimation method.
sigmaValues_TiGo <- rep(NA,161)
n_drop <- c()
par_vec <- matrix(NA,161,4)
# Set alpha and beta
alpha <- 0.2
beta <- 0.04
for (j in 1:161) {
  if (is.na(startend_dates$end_year_ts[j])) next
  y1 <- as.numeric(timeSeriesSearchInterests[,j])
  y <- c(na.omit(y1))
  y[y<=0] <- min(y[y>0])/2
  TiGoModel <- TiGo.ETS(y,alpha=alpha,beta=beta)
  par_vec[j,] <- TiGoModel$par[3:6] # store parameters
  sigmaValues_TiGo[j] <- TiGoModel$sigma # Store the sd of the errors to sigma2Values
  # If the model suggests that the product has not passed its peak, we remove this product from the ensemble
  if (tail(exp(TiGoModel$states[,2]),1)>1) n_drop <- c(n_drop,j)
}

# Set up the accuracy matrix
accuracy_result <- array(NA,dim = c(161,209,24,99))
n_drop <- unique(n_drop)
TiGoModel1 <- TiGoModel

for (i in 1:161) {
  y1 <- as.numeric(timeSeriesSearchInterests[,i])
  y <- na.omit(y1)
  H1 <- length(y) - 12
  H2 <- length(y) - 24
  if (H1 <= 0) next
  for (j in 0:H1) {
    n_prior <- which((end < start[i]+months(j-1)))
    n_prior <- n_prior[!n_prior%in%n_drop]
    if (length(n_prior) < 2) next
    ytrain <- y[1:j]
    if (j <= H2)  ntest <- 24
    if (j > H2 & j <= H1)  ntest <- 12
    # We replace zeros in these series with the minimum positive value in the series, divided by 2.
    # If all the data points up to t are zeros, we replace them with the minimum positive value in the completed life cycles, divided by 2.
    if (any(ytrain<=0)) {
      if (all(ytrain == 0)) {
        prior_curves <- timeSeriesSearchInterests[n_prior,]
        prior_curves[prior_curves == 0] <- 100
        ytrain[which(ytrain==0)] <- min(prior_curves,na.rm=TRUE)/2
      } else {
        ytrain[which(ytrain==0)] <- min(ytrain[ytrain>0])/2
      }
    }
    ytrain <- log(ytrain)
    
    # Use the mean and variance of the parameters estimates to form prior distributions.
    phi <- mean(par_vec[n_prior,1])
    tau <- mean(par_vec[n_prior,4])
    l0 <- mean(par_vec[n_prior,2])
    b0 <- mean(par_vec[n_prior,3])
    sigmam_TiGo <- mean(sigmaValues_TiGo[n_prior], na.rm = TRUE)
    if (length(n_prior) > 2) {
      phiValues1 <- par_vec[n_prior,1][which(par_vec[n_prior,1] <= quantile(par_vec[n_prior,1],0.75)&par_vec[n_prior,1]>=quantile(par_vec[n_prior,1],0.25))]
      tauValues1 <- par_vec[n_prior,4][which(par_vec[n_prior,4] <= quantile(par_vec[n_prior,4],0.75)&par_vec[n_prior,4]>=quantile(par_vec[n_prior,4],0.25))]
      l0Values1 <- par_vec[n_prior,2][which(par_vec[n_prior,2] <= quantile(par_vec[n_prior,2],0.75)&par_vec[n_prior,2]>=quantile(par_vec[n_prior,2],0.25))]
      b0Values1 <- par_vec[n_prior,3][which(par_vec[n_prior,3] <= quantile(par_vec[n_prior,3],0.75)&par_vec[n_prior,3]>=quantile(par_vec[n_prior,3],0.25))]
    } else {
      phiValues1 <- par_vec[n_prior,1]
      tauValues1 <- par_vec[n_prior,4]
      l0Values1 <- par_vec[n_prior,2]
      b0Values1 <- par_vec[n_prior,3]
    }
    if (length(phiValues1) == 1 | length(tauValues1) == 1 | length(b0Values1) == 1 | length(l0Values1) == 1) {
      phiValues1 <- par_vec[n_prior,1]
      tauValues1 <- par_vec[n_prior,4]
      l0Values1 <- par_vec[n_prior,2]
      b0Values1 <- par_vec[n_prior,3]
    }
    initial.cov.matrix_TiGo <- cov(cbind(na.omit(l0Values1),na.omit(b0Values1),na.omit(phiValues1),na.omit(tauValues1)))
    l0v<-initial.cov.matrix_TiGo[1,1]
    b0v<-initial.cov.matrix_TiGo[2,2]
    phiv<-initial.cov.matrix_TiGo[3,3]
    tauv<-initial.cov.matrix_TiGo[4,4]
    
    # Set up variables and matrices for the Kalman filter
    F <- matrix(c(1, 0, phi, phi), 2, 2)
    w <- matrix(c(1, phi), 2, 1)
    A <- rbind(t(w), F)
    g <- matrix(c(alpha, beta), 2, 1)
    b <- matrix(c(1, g), 3, 1)
    c31 <- t(t(c(0,0,log(tau))))
    m_t1t1 <- matrix(c(l0, b0), 2, 1)
    V_t1t1 <- matrix(c(abs(l0v), 0,0, abs(b0v)), 2, 2) 
    
    oneStep <- matrix(NA, j,1)
    z <- matrix(NA, j, 2)
    
    if (j >= 1) {
      for (t in 1:j) {
        # Prediction step
        mz_tt1 <- A%*%m_t1t1 + c31
        mu <- mz_tt1[1]
        Vz_tt1 <- A%*%V_t1t1%*%t(A) + sigmam_TiGo^2*b%*%t(b)
        v_tt1 <- Vz_tt1[1,1]
        V_tt1 <- Vz_tt1[2:3, 2:3]
        eta <- Vz_tt1[2:3, 1]
        
        # Revision step
        k_t <- eta/v_tt1
        V_t1t1 <- V_tt1 - v_tt1*k_t%*%t(k_t)
        m_t1t1 <- F%*%m_t1t1 + k_t*(ytrain[t] - mu)
        
        oneStep[t] <- mu
        z[t, ] <- m_t1t1
      }
    } else {
      z <- c(l0,b0)
    }
    
    # Generate quantile forecasts from the updated parameters
    par <- c(alpha,beta,phi,l0,b0,tau)
    TiGoModel1$par <- par
    TiGoModel1$y <- exp(ytrain)
    if (length(z) == 2) {
      TiGoModel1$states <- matrix(z,1,2)
    } else {
      TiGoModel1$states <- z
    }
    
    TiGoModel1$sigma <- sigmam_TiGo
    ytest <- y[(j+1):(j+ntest)]
    forecastTp <- forecast.TiGo.ETS(TiGoModel1,h=ntest,level=seq(2,98,2))
    forecastTp$predictions <- forecastTp$predictions[,-c(2,3)]
    forecastTp$predictions[forecastTp$predictions < 0] <- 0 
    forecastTp$predictions[forecastTp$predictions > 100] <- 100 
    accuracy_result[i,j+1,1:ntest, 50] <-  pinball(forecastTp$predictions[1:ntest,1], ytest,0.5)
    
    # Calcuate pinball losses
    for (mm in 1:99) {
      if (mm <= 49) {
        accuracy_result[i,j+1,1:ntest,mm] <- pinball(forecastTp$predictions[1:ntest,100-2*mm], ytest,quant[mm])
      }
      if (mm >= 51) {
        accuracy_result[i,j+1,1:ntest,mm] <- pinball(forecastTp$predictions[1:ntest,2*mm-99], ytest,quant[mm])
      }
    }
  }
}

# Print out average accuracy results
result_vec <- t(as.matrix(2*c(mean(rowMeans(apply(accuracy_result[,,1:12,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                              mean(rowMeans(apply(apply(accuracy_result[,,1:12,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                              mean(rowMeans(apply(accuracy_result[,,13:24,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                              mean(rowMeans(apply(apply(accuracy_result[,,13:24,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE))))
colnames(result_vec) <- c('1-12 steps MAE','1-12 steps CRPS','13-24 steps MAE','13-24 steps CRPS')
result_vec
