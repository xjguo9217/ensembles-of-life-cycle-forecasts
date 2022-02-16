# -----------------------------------------------------------------
# Forecasting Search Interests of Social Networking Websites
# Using Machine Learning Models
# -----------------------------------------------------------------

# In this script, we use Quantile Regression Forests and LightGBM to predict search interests of social networking websites

# Choose a model
model_select <- "QRF" # This can be set to QRF or LightGBM

# Load in the lubridate package
library(lubridate)
options(warn=-1)

# Load in the data of time series and start/end dates of each time series
data1 <- read.csv('GoogleTrendsTimeSeriesSocial-May-2021.csv')
startend_dates <- read.csv('WikiListSocial-startenddays-May-2021.csv') 

timeSeriesSearchInterests <- data1[,-c(1:2)] # remove the time columns
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

# Find unique sets of completed life cycles.
priors_list <- c()
for (i in 1:161) {
  y1 <- as.numeric(timeSeriesSearchInterests[, i])
  y <- na.omit(y1)
  H1 <- length(y) - 12
  if (H1 <= 0) next
  for (j in 0:H1) {
    n_complete <- which((end < start[i]+months(j-1)))
    if (length(n_complete) < 2) next
    priors_list <- c(priors_list,list(which((end < start[i]+months(j-1)))))
  }
}
priors_list_unique <- unique(priors_list)
priors_list <- c()

# Train QRF/LightGBM model using times series in each set of completed life cycles.
model_list <- c()
window <- 24 # We use y_{t-24,t-1} as predictors to predict y_t 
for (j in 1:length(priors_list_unique)) {
  training_data <- rep(NA,window + 1)
  n1 <- priors_list_unique[[j]]
  for (i in n1) {
    y1 <- as.numeric(timeSeriesSearchInterests[,i])
    y <- as.numeric(na.omit(y1))  
    for (kk in 1:length(y)) {
      if (kk == 1) this_t <-c(rep(0,window),y[1])
      if (kk <= window & kk > 1) this_t <- c(rep(0,window+1-kk),y[1:kk]) 
      if (kk > window) this_t <- y[(kk-window):kk]
      training_data <- rbind(training_data,this_t)
    }
  }
  training_data <- training_data[-1,]
  set.seed(201)
  if (model_select == 'QRF') {
    qrf  <- quantregForest(training_data[,1:window], training_data[,window+1])
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
      dtrain <- lgb.Dataset(training_data[,1:window], label = training_data[,window+1])
      model_list <- c(model_list,lgb.train(params = lgb.grid, data = dtrain))
    } 
  }
}

# Set up the accuracy matrix
accuracy_result <- array(NA,dim = c(161,202,24,99))

# Generate predictions
for (i in 1:161) {
  y1 <- as.numeric(timeSeriesSearchInterests[,i])
  y <- na.omit(y1)
  H1 <- length(y) - 12
  H2 <- length(y) - 24
  if (H1 < 0) next
  # Generating rolling forecasts from t=0 to t=H1
  for (j in 0:H1) {
    # Find products that are completed before t
    ytrain <- y[1:j]
    n_complete <- which((end < start[i]+months(j-1)))
    if (length(n_complete) < 2) next # Move to the next one if there are less than two completed life cycles
    # Find the model corresponding to the set of completed life cycles
    n_model <-  which(lapply(priors_list_unique,function(x) length(x) == length(n_complete) && all(x == n_complete)) == TRUE)
    qrf_model <- model_list[n_model]
    if (j <= H2)  ntest <- 24
    if (j > H2 & j <= H1)  ntest <- 12
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
      if (model_select == "QRF") forecastTp[mm,] <-  predict(qrf_model[[1]], this_test, what = seq(0.01,0.99,0.01)) 
      if (model_select == "LightGBM") {
        for (cc in 1:99) {
          gbm_model <- model_list[(n_model-1)+cc]
          this_test <- rbind(this_test)
          forecastTp[mm,cc] <-  predict(gbm_model[[1]], this_test) 
        }
      }
    }
    forecastTp[forecastTp>100] <- 100
    
    # Calculate pinball losses
    for (mm in 1:99) {
      accuracy_result[i,j+1,1:ntest,mm] <- pinball(forecastTp[,mm], ytest, mm/100)
    }
  }
}

# Print out average accuracy results
result_vec <- t(as.matrix(2*c(mean(rowMeans(apply(accuracy_result[,,1:12,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                  mean(rowMeans(apply(apply(accuracy_result[,,1:12,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                  mean(rowMeans(apply(accuracy_result[,,13:24,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                  mean(rowMeans(apply(apply(accuracy_result[,,13:24,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE))))
colnames(result_vec) <- c('1-8 steps MAE','1-8 steps CRPS','9-17 steps MAE','9-17 steps CRPS')
result_vec
