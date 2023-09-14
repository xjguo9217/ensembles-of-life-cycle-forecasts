# --------------------------------------------------------------
# Forecasting Demand of Intel Microprocessors
# Using Time Varying Bass Diffusion Model Updated Using Extended 
# Kalman Filter with Continuous States and Discrete Observations
# ---------------------------------------------------------------

# In this script, we generate forecasts of Intel microprossors sales using the Bass-Extended Kalman filter with continuous state and discrete observations method.

# Load in the lubridate and deSolve package
library(deSolve)
library(zoo)
source('Source code/Bass.R') # Load in functions for the Bass model.
options(warn=-1)

# Load in the data of time series and start/end dates of each time series
data1 <- read.csv('Intel_actual_demand.csv')
startend_date <- read.csv('startend_date.csv')

timeSeriesIntel <- data1[,-c(1)] # remove the time columns
quant <- quant <- seq(0.01,0.99,0.01) # we consider 99 quantiles in model evaluation.

# Set up a function for calculating the pinball loss
pinball <- function(pred, act, q) {
  I = 1*(act>pred)
  l <- (act - pred)*q*I + (pred - act)*(1-q)*(1-I)
  return(l)
}

start <- startend_date$start
end <- startend_date$end

# Find unique sets of completed life cycles.
priors_list <- c()
for (i in 1:86) {
  y1 <- as.numeric(timeSeriesIntel[, i])
  y <- na.omit(y1)
  H1 <- length(y) - 8
  if (H1 <= 0) next
  for (j in 0:H1) {
    n_complete <- which(end < start[i]+j-1)
    if (length(n_complete) < 2) next
    priors_list <- c(priors_list,list(which((end < start[i]+j-1))))
  }
}
priors_list_unique <- unique(priors_list)
priors_list <- c()

# Fit the Bass model to all time series using the Maximum Likelihood Estimation method.
# Use the median and variance of the parameters estimates to form prior distributions.
M_values <- rep(NA,length(priors_list_unique))
p_values <- rep(NA,length(priors_list_unique))
q_values <- rep(NA,length(priors_list_unique))
var_values <- rep(NA,length(priors_list_unique))
pv_values<-rep(NA,length(priors_list_unique))
qv_values<-rep(NA,length(priors_list_unique))
mv_values<-rep(NA,length(priors_list_unique))
n_drop_all <- c()
for (j in 1:length(priors_list_unique)) {
  n1 <- priors_list_unique[[j]]
  n_drop <- c()
  sigma2Values_bass <- rep(NA,length(n1))
  n_length <- c()
  for (i in n1) {
    y1 <- as.numeric(timeSeriesIntel[,i])
    y <- as.numeric(na.omit(y1))  
    n_length <- c(n_length,length(y))
  }
  n_length <- max(n_length)
  fitCurves <- matrix(NA,length(n1), n_length)
  sigmacum2Values_bass <- rep(NA, length(n1))
  mValues <- rep(NA, length(n1))
  pValues <- rep(NA, length(n1))
  qValues <- rep(NA, length(n1))
  for (i in 1:length(n1)) {
    y1 <- as.numeric(timeSeriesIntel[,n1[i]])
    y <- as.numeric(na.omit(y1))  
    BassModel <- Bass.fit(y,error='normal')
    sigma2Values_bass[i] <- BassModel$sigma^2
    if (which.max(BassModel$fit.y) == length(BassModel$fit.y)) {
      n_drop <- c(n_drop, i)
      n_drop_all <- c(n_drop_all, n1[i])
    }
    BassModel$y <- NULL
    fitCurves[i,1:n_length] <- forecast.Bass(BassModel,n_length,error="normal")$prediction[,1]
    mValues[i] <- BassModel$par[3]
    pValues[i] <- BassModel$par[1]
    qValues[i] <- BassModel$par[2]
  }
  if (length(n_drop)>0) {
    fitCurves <- fitCurves[-n_drop,]
    mValues <- mValues[-n_drop]
    pValues <- pValues[-n_drop]
    qValues <- qValues[-n_drop]
    sigma2Values_bass <- sigma2Values_bass[-n_drop]
  }
  M_values[j] <- mean(mValues,na.rm=TRUE, trim = 0.5)
  p_values[j] <- mean(pValues,na.rm=TRUE, trim = 0.5)
  q_values[j] <- mean(qValues,na.rm=TRUE, trim = 0.5)
  var_values[j] <- mean(sigma2Values_bass, na.rm = TRUE, trim = 0.5)
  # We use only values within the interquartile range to calculate variances.
  if (length(pValues) > 2) {
    pValues <- pValues[which(pValues <= quantile(pValues,0.75)&pValues>=quantile(pValues,0.25))]
    qValues <- qValues[which(qValues <= quantile(qValues,0.75)&qValues>=quantile(qValues,0.25))]
    mValues <- mValues[which(mValues <= quantile(mValues,0.75)&mValues>=quantile(mValues,0.25))]
  }
  initial.cov.matrix_bass <- cov(cbind(na.omit(pValues),na.omit(qValues),na.omit(mValues)))
  pv_values[j]<-initial.cov.matrix_bass[1,1]
  qv_values[j]<-initial.cov.matrix_bass[2,2]
  mv_values[j]<-initial.cov.matrix_bass[3,3]
}
n_drop_all = unique(n_drop_all)

# Set up accuracy matrix
accuracy_result <- array(NA,dim = c(86,187,17,99))
set.seed(201)
for (i in 1:86) {
  y1 <- as.numeric(timeSeriesIntel[,i])
  y <- na.omit(y1)
  H1 <- length(y) - 8
  H2 <- length(y) - 17
  if (H1 < 0) next
  n_change <- c()
  priors_list <- c(0,0,0,0,0)
  # For time series i, find all the time points when prior information fold changes
  # These changing points will divided the time series into several periods
  for (j in 0:(H1)) {
    n_prior <- which(end < start[i]+j-1)
    if (length(n_prior) < 2) next
    if (!all(n_prior%in%priors_list)) {
      priors_list<-n_prior
      n_change <- c(n_change,j)
    }
  }
  # For each period between two changing points, generate forecasts from the Bass-EKF model
  for (aaa in 2:(length(n_change)+1)) {
    if (aaa == (length(n_change)+1)) {
      nn <- H1
    } else {
      nn <- (n_change[aaa]-1)
    }
    n_prior <- which((end < start[i]+nn-1))
    if (length(n_prior) < 2) next
    n_model <-  which(lapply(priors_list_unique,function(x) length(x) == length(n_prior) && all(x == n_prior)) == TRUE)
    # Find the mean and variance of the parameters
    p <- p_values[n_model]
    q <- q_values[n_model]
    m <- M_values[n_model]
    sigma2m_bass <- var_values[n_model]
    pv_bass <- pv_values[n_model]
    qv_bass <- qv_values[n_model]
    mv_bass <- mv_values[n_model]
    
    start_t <- 1
    srat_N <- 0
    # Set up a function to calculate dN/dt and dP/dt in the differential equations
    modelFunction<- function(Time, State, Pars) {
      with(as.list(c(State, Pars)), {
        a <- q-p-2*q*N/m
        b <- m - N
        c <- N - N^2/m
        d <- p + q*N^2/m^2
        dN <- (p+q/m*N)*(m-N)
        dP1 <- 2*a*P1 + 2*b*P2 + 2*c*P3 + 2*d*P4 + sigma2
        dP2 <- a*P2 + b*u + c*h + d*jj
        dP3 <- a*P3+ b*h + c*v + d*k
        dP4 <- d*w + a*P4 + b*jj + c*k
        return(list(c(dN, dP1, dP2, dP3, dP4)))
      })
    }
    
    # Set up a function to calculate dN/dt
    modelFunction1 <- function(Time, State, Pars) {
      with(as.list(c(State, Pars)), {
        dN <- (p+q/m*N)*(m-N)
        return(list(c(dN)))
      })
    }
    
    # Initialize the filter at t=0
    priorBeta <- c(p, q, m)
    H <- c(1, 0, 0, 0)
    w1 <- sigma2m_bass   
    Q <-  matrix(c(w1, rep(0,15)), 4, 4) 
    s0 <- 0
    P0 <- matrix(c(s0, 0, 0, 0, 0, pv_bass, 0, 0, 0, 0, qv_bass, 0, 0,0,0,mv_bass), 4, 4)
    I <- diag(4)
    y_post <- t(t(c(s0,priorBeta)))
    P_post <- P0
    
    # Run Extended Kalman Filter to generate forecasts
    fit  <- c()
    K <- sigma2m_bass
    for (t in 1:(nn+1)) {
      if (t <= H2+1)  ntest <- 17
      if (t > H2+1 & t <= H1+1)  ntest <- 8
      ytest <- y[(t):(t+ntest-1)]
      step <- ntest
      # Initialize the filter
      pt <- y_post[2]; qt <- y_post[3]; mt <- y_post[4]
      n <- y_post[1]
      
      # Predictive steps
      # Solve differential equations to update sales forecast N and covariance of the state P
      beta_post <- matrix(rep(NA,(nn + 1)*3),(nn  + 1),3)
      speriod <- rep(NA, nn) 
      fitsumy <- c()
      pars <- c(p = pt, q = qt, m = mt, u = P_post[2,2],
                v = P_post[3,3], w = P_post[4,4], h = P_post[2,3], 
                jj = P_post[2,4], k = P_post[3,4],sigma2 = Q[1,1])
      yini <- c(N = n, P1 = P_post[1,1], P2 = P_post[1,2], P3 = P_post[1,3], P4 = P_post[1,4])
      times <- seq(c(t-1+start_t), c(t+start_t), by = 0.001)
      ntime <- 1/0.001+1
      out <- ode(func = modelFunction, y = yini, parms = pars, times = times,method='rk')
      if (out[ntime,2] == Inf) ntime <- 1
      
      # Collect parameters for generating quantile forecasts
      s <- out[ntime,2]
      P1 <- out[ntime,3]
      P2 <- out[ntime,4]
      P3 <- out[ntime,5]
      P4 <- out[ntime,6]
      P_prior <- P_post
      P_prior[1,] <- c(P1,P2,P3,P4)
      P_prior[,1] <- c(P1,P2,P3,P4)
      speriod[t] <- s - n 
      y_prior <- t(t(c(s, y_post[2:4])))
      npred <- c(s,rep(NA,step))
      speriodpred <- rep(NA,step)
      speriodpred[1] <- speriod[t]
      fit  <- c(fit,speriodpred[1])
      
      # Generate quantile forecasts
      if (t-1 >= n_change[aaa-1]) {
        speriodpred_q <- matrix(NA,step,99)
        for (a in 1:99) {
          speriodpred_q[1,a] <- speriod[t] + qnorm(a/100)*sqrt(2)*sqrt(K)
        }
        for (j in 2:step) {
          pars <- c(p = pt, q = qt, m = mt)
          yini <- c(N = as.numeric(npred[j-1]))
          times <- seq(c(t+start_t+j-2), c(t+start_t+(j-1)), by = 0.001)
          ntime <- 1/0.001+1
          out <- ode(func = modelFunction1, y = yini, parms = pars, times = times,method='rk')
          if (out[ntime,2] == Inf) ntime <- 1
          npred[j] <- out[ntime,2]
          speriodpred[j] <- max(0,npred[j] - npred[j-1])
          for (a in 1:99) {
            speriodpred_q[j,a] <- speriodpred[j] + qnorm(a/100)*sqrt(2)*sqrt(K)
          }
        }
        speriodpred_q[speriodpred_q<0] <- 0 # replace forecasts < 0 with 0
        
        # Calculate pinball losses
        for (mm in 1:99) {
          accuracy_result[i,t,1:ntest,mm] <- pinball(speriodpred_q[,mm], ytest, mm/100)
        }
      }
      # Update steps
      R <- sigma2m_bass
      phi <- P_prior%*%t(t(H))%*%(1/(H%*%P_prior%*%t(t(H)) + R)) # calculate Kalman Gain
      if (!any(is.na(phi))) {
        y_post <- y_prior + phi%*%(sum(y[1:(t)]) - H%*%y_prior)   
        P_post <- (I-phi%*%H)%*%P_prior
        K <- H%*%P_prior%*%t(t(H)) + R
      } 
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
