# ---------------------------------------------
# Forecasting Demand of Intel Microprocessors
# Using Machine Learning Models
# ---------------------------------------------

# In this script, we use Quantile Regression Forests and LightGBM to predict sIntel microprocessors sales.
library(Metrics)
# Choose a model
model_select <- "QRF" # This can be set to QRF or LightGBM

options(warn=-1)

# Load in the data of time series and start/end dates of each time series
data1 <- read.csv('Intel_actual_demand.csv')
startend_date <- read.csv('startend_date.csv')

timeSeriesIntel <- data1[,-c(1)] # remove the time columns
quant <- quant <- seq(0.01,0.99,0.01) # we consider 99 quantiles in model evaluation.

# Load in functions for the selected candidate model.
if (model_select == 'QRF') library(quantregForest)
if (model_select == 'LightGBM') library(lightgbm) 

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

# Train QRF/LightGBM model using times series in each set of completed life cycles.
model_list <- c()
window <- 17 # We use y_{t-17,t-1} as predictors to predict y_t 
for (j in 1:length(priors_list_unique)) {
  training_data <- rep(NA,window + 1)
  n1 <- priors_list_unique[[j]]
  for (i in n1) {
    y1 <- as.numeric(timeSeriesIntel[,i])
    y <- as.numeric(na.omit(y1))  
    for (kk in 1:length(y)) {
      if (kk == 1) this_t <-c(rep(0,window),y[1])
      if (kk <= window & kk > 1) this_t <- c(rep(0,window+1-kk),y[1:kk]) 
      if (kk > window) this_t <- y[(kk-window):kk]
       training_data <- rbind(training_data,this_t)
    }
  }
  training_data <- data.frame(training_data[-1,])
  
  set.seed(201)
  if (model_select == 'QRF') {
    qrf  <- quantregForest(training_data[,-(window+1)], as.numeric(training_data[,window+1]),ntree=250)
    model_list <- c(model_list,list(qrf))
  }
  if (model_select == 'LightGBM') {
    for (cc in 1:99) {
      lgb.grid = list(objective = "quantile",
                      metric = "quantile",
                      alpha=quant[cc],
                      learning_rate= 0.075,
                      sub_row= 0.75,
                      bagging_freq= 1,
                      lambda_l2= 0.1)
      dtrain <- lgb.Dataset(as.matrix(training_data[,-(window+1)]), label = as.numeric(training_data[,window+1]))
      model_list <- c(model_list,lgb.train(params = lgb.grid, data = dtrain))
    } 
  }
}

# Set up the accuracy matrix
accuracy_result <- array(NA,dim = c(86,187,17,99))
# Generate predictions
par(mfrow=c(3,4))
for (i in 1:86) {
  y1 <- as.numeric(timeSeriesIntel[,i])
  y <- na.omit(y1)
  H1 <- length(y) - 8
  H2 <- length(y) - 17
  if (H1 < 0) next
  # Generating rolling forecasts from t=0 to t=H1
  for (j in 0:H1) {
    # Find products that are completed before t
    ytrain <- y[1:j]
    n_complete <- which((end < start[i]+j-1))
    if (length(n_complete) < 2) next # Move to the next one if there are less than two completed life cycles
    # Find the model corresponding to the set of completed life cycles
    n_model <-  which(lapply(priors_list_unique,function(x) length(x) == length(n_complete) && all(x == n_complete)) == TRUE)
    qrf_model <- model_list[n_model]
    if (j <= H2)  ntest <- 17
    if (j > H2 & j <= H1)  ntest <- 8
    ytest <- y[(j+1):(j+ntest)]
    # Create predictors in the test set
    if (j  == 0) this_test <-c(rep(0,window))
    if (j <= window & j  > 0) this_test <- c(rep(0,window-(j)),y[1:(j)]) 
    if (j > window) this_test <- y[c(j-window+1):(j)]
    forecastTp <- matrix(NA,ntest,99)
    for (mm in 1:ntest) {
      if (mm == 1) {
         this_test <- this_test
      } else {
         this_test <- c(this_test[2:window],forecastTp[mm-1,50])
      }
  #    if (i == 150 & j==54) this_test <- this_test + 2
      if (model_select == "QRF") forecastTp[mm,] <-  predict(qrf_model[[1]], this_test, what = seq(0.01,0.99,0.01)) 
      if (model_select == "LightGBM") {
        for (cc in 1:99) {
          gbm_model <- model_list[(n_model-1)+cc]
          this_test <- rbind(this_test)
          forecastTp[mm,cc] <-  predict(gbm_model[[1]], this_test) 
        }
      }
    }
    
    # Calculate pinball losses
    for (mm in 1:99) {
      accuracy_result[i,j+1,1:ntest,mm] <- pinball(forecastTp[,mm], ytest, mm/100)
    }
  }
}

# split to train/test data
start_allts <- 1
end_allts <- 187
n_train <- round(187/2)
end_train <- start_allts + n_train
end_train_ts <- end_train - start +1

accuracy_train <- array(NA,dim = c(86,187,17,99))
accuracy_test <- array(NA,dim = c(86,187,17,99))
for (i in 1:86) {
  this_end_train <- end_train_ts[i]
  if (this_end_train <= 1) {
    accuracy_test[i,,,] <- accuracy_result[i,,,] 
  } else {
    accuracy_train[i,1:(this_end_train-1),,] <- accuracy_result[i,1:(this_end_train-1),,] 
    accuracy_test[i,this_end_train:187,,] <- accuracy_result[i,this_end_train:187,,] 
  }
}

result_test <- t(2*c(mean(rowMeans(apply(accuracy_test[,,1:8,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                     mean(rowMeans(apply(apply(accuracy_test[,,1:8,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                     mean(rowMeans(apply(accuracy_test[,,9:17,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                     mean(rowMeans(apply(apply(accuracy_test[,,9:17,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE)))
colnames(result_test) <- c("1-8 steps MAE","1-8 steps MCRPS","9-17 steps MAE","9-17 steps MCRPS")

result_train <- t(2*c(mean(rowMeans(apply(accuracy_train[,,1:8,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                      mean(rowMeans(apply(apply(accuracy_train[,,1:8,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                      mean(rowMeans(apply(accuracy_train[,,9:17,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                      mean(rowMeans(apply(apply(accuracy_train[,,9:17,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE)))
colnames(result_train) <- c("1-8 steps MAE","1-8 steps MCRPS","9-17 steps MAE","9-17 steps MCRPS")

result_vec
result_test
