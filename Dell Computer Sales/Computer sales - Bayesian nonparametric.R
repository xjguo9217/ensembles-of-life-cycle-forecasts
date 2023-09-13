# -------------------------------------------------------------------------
# Forecasting Computer Sales Using the Bayesian Ensemble of TiGo-ETS Models
# Using the Bayesian non-parametric method 
# -------------------------------------------------------------------------
# The Stan Code is adapted from the script provided in Dew et al.'s (2018) Web Appendix C 

# Load in the package
library(boot)
library(rstan)
library(foreach)
library(doSNOW)
library(doParallel)
library(Metrics)
library(lubridate)
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
CLUSTER_SIZE=3
cl <- makeCluster(3)
registerDoParallel(cl)

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

accuracy_result <- array(NA,dim = c(170,82,34,99))

# Stan code from Web Appendix C of Dew et al. (2018)
stanmodelcode <- '
data{
  int<lower=1> P; // number of periods
  int<lower=1> pred_h; // number of periods
  int<lower=1> N; // number of customers
  int<lower=1> M; // number of data points
  int<lower=1> id[M]; // customer id of observation m
  int<lower=1> t[M]; // calendar time
  int<lower=0,upper=1> y[M-pred_h]; // outcome of observation m (purchase or not)
}
parameters{
  // non-restricted components of the GP
vector[P] alpha_long;
// unobs. heterogeneity, random effects
vector[N] delta;
real<lower=0> sigsq;
// kernel hyperparameters
real<lower=0> etasq_long;
positive_ordered[2] rhosq_cal;

// mean function hyperparameters
real mu; // baseline spending level (cal mean function)
}
transformed parameters{
// full GPs (with the first period restriction)
vector[P] long_mean;

vector[1] z;
// add zero first period restriction
z[1] <- 0;

// set mean functions
for(p in 1:P) {
long_mean[p] <- mu;}
}
model{
vector[M] theta;
matrix[P,P] Sigma_long;
matrix[P,P] L_long;

// calendar time, long-run component kernel
for(i in 1:P) {
for(j in 1:P) {
Sigma_long[i,j] <- etasq_long * exp(-((i-j)^2)/rhosq_cal[2]);
Sigma_long[j,i] <- Sigma_long[i,j]; }
Sigma_long[i,i] <- Sigma_long[i,i]+0.0001; }
L_long <- cholesky_decompose(Sigma_long);

// kernel hyperparameter priors
etasq_long ~ normal(0,5);
rhosq_cal ~ normal(P/2,P);

// kernel mean function priors
mu ~ normal(0,5);

// sample GPs
alpha_long ~ multi_normal_cholesky(long_mean,L_long);

// sample random effects
delta ~ normal(0,sigsq);
sigsq ~ normal(0,2.5);
// gppm likelihood
for(m in 1:M)
theta[m] <- alpha_long[t[m]]+delta[id[m]];
y ~ bernoulli_logit(theta[1:(M-pred_h)]);
}
'
stanModel <- stan_model(model_code = stanmodelcode)

# Produce forecasts
for (i in 1:170) {
  y1 <- as.numeric(timeSeriesDell[i,])
  y <- na.omit(y1)
  H1 <- length(y) - 8
  H2 <- length(y) - 17
  if (H1 <= 0) next
  for (j in 0:H1){
    ntrain <- j
    n_prior <- which((end < start[i]+days(7*(j-1))))
    if (length(n_prior) < 2) next
    ytrain <- y[1:j]
    if (j <= H2)  ntest <- pred_h <- 17
    if (j > H2 & j <= H1)  ntest <- pred_h <- 8
    ytest <- y[(ntrain + 1):(ntrain + ntest)]
    maxm <- c()
    for (l in 1:length(n_prior)) {
      y11 <- as.numeric(timeSeriesDell[n_prior[l],])
      maxm <- c(maxm,sum(na.omit(y11)))
    }
    # Assume the number of potential customers is max(the median total number of customers for completed life cycles, sum of the current total number of customers in this series).
    mm = max(round(median(maxm)),sum(ytrain)) # scale down the number of potential customers to make the algorithm faster.
    mm1 <- round(mm/(max(ytrain)/10))
    if (mm1 == Inf) mm1 <- 100
    # Set up the Bayesian model.
    yy = rep(0,length(ytrain)*mm1)
    ytrain1 <- ytrain/mm*mm1
    for (kk in 1:length(ytrain)) {
      if(1 > (round(ytrain1[kk]))) next
      yy[min((kk-1)*mm1 + length(which(yy == 1))+1,kk*mm1) : min(((kk-1)*mm1 + length(which(yy == 1))+round(ytrain1[kk])),kk*mm1)] <- 1
    }
    tt <- rep(0,(length(ytrain)+pred_h)*mm1)
    for (kk in 1:(length(ytrain)+pred_h)) {
      tt[((kk-1)*mm1+1) : (kk*mm1)] <- kk
    }
    id <- rep(1:mm1,(length(ytrain)+pred_h))
    data <- list(
      N = mm1,
      M = (length(ytrain)+pred_h)*mm1,
      t = tt,
      id = id,
      pred_h = pred_h*mm1,
      P = (length(ytrain)+pred_h),
      y = yy
    )
    ret<- foreach(i=1:CLUSTER_SIZE, .packages=c("rstan")) %dopar% {#run sampling in parallel
      sampling(stanModel,
               data=data,
               iter=100,  # increase the number of iterations to get better results.
               chains=2,
               refresh = 0,verbose=TRUE)
    }
    samples <- sflist2stanfit(ret)
    inv_log <- extract(samples)$alpha_long#+extract(samples)$alpha_short #extract(samples)$alpha_week++ extract(samples)$delta
    if (j > 1) {
      fit <- inv.logit(apply(inv_log[,1:ntrain],2,median)[1:ntrain])*mm
    } else {
      fit <- inv.logit(median(inv_log[,1:ntrain])[1:ntrain])*mm
    }
    samples <- c()
    ret <- c()
    inv_log <- inv_log[,(ntrain+1):(ntrain+ntest)]
    prob <- t(inv.logit(apply(inv_log,2,quantile,probs=quant)))
    quantile_pred <- prob*mm
  
    # calculate accuracy
    for (mm in 1:99) {
      accuracy_result[i,j+1, 1:ntest, mm] <- pinball(quantile_pred[,mm], ytest,quant[mm])
    }
  }  
}

# Print out average accuracy results
accuracy_result[which(accuracy_result == Inf)] <- NA
result_vec <- t(2*c(mean(rowMeans(apply(accuracy_result[,,1:8,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                    mean(rowMeans(apply(apply(accuracy_result[,,1:8,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                    mean(rowMeans(apply(accuracy_result[,,9:17,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                    mean(rowMeans(apply(apply(accuracy_result[,,9:17,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE)))
colnames(result_vec) <- c("1-8 steps MAE","1-8 steps MCRPS","9-17 steps MAE","9-17 steps MCRPS")
result_vec
