# -----------------------------------------------------------------
# Forecasting Computer Sales Using Machine Learning Models
# -----------------------------------------------------------------

# In this script, we use Quantile Regression Forests and LightGBM to predict Dell computer sales.

# Choose a model
model_select <- "QRF" # This can be set to QRF or LightGBM

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

# Load in functions for the selected candidate model.
if (model_select == 'QRF') library(quantregForest)
if (model_select == 'LightGBM') library(lightgbm) 

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

# Find unique sets of completed life cycles.
priors_list <- c()
for (i in 1:170) {
  y1 <- as.numeric(timeSeriesDell[i, ])
  y <- na.omit(y1)
  H1 <- length(y) - 7
  if (H1 <= 0) next
  for (j in 1:H1) {
    n_prior <- which((end < start[i]+days(7*(j-1))))
    if (length(n_prior) < 2) next
    priors_list <- c(priors_list,list(which((end < start[i]+days(7*(j-1))))))
  }
}
priors_list_unique <- unique(priors_list)
priors_list <- c()

# Train the QRF/LightGBM model using times series in each set of completed life cycles.
model_list <- c()
window <- 34 # We use y_{t-35,t-1} as predictors to predict y_t 
for (j in 1:length(priors_list_unique)) {
  training_data <- rep(NA,window + 1)
  n1 <- priors_list_unique[[j]]
  for (i in n1) {
    y1 <- as.numeric(timeSeriesDell[i, ])
    y <- as.numeric(na.omit(y1))  
    for (j in 1:length(y)) {
      if (j == 1) this_t <-c(rep(0,window),y[1])
      #  if j < 34, the time series is padded upfront by 35-j zeros for pre-launch forecasting.
      if (j <= window & j > 1) this_t <- c(rep(0,window+1-j),y[1:j]) 
      if (j > window) this_t <- y[(j-window):j]
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
accuracy_result <- array(NA,dim = c(170,82,34,99))

# Generate predictions
for (i in 1:170) {
  y1 <- as.numeric(timeSeriesDell[i, ])
  y <- na.omit(y1)
  H1 <- length(y) - 8
  H2 <- length(y) - 17
  if (H1 < 0) next
  # Generate rolling forecasts from t=0 to t=H1
  for (j in 0:H1) {
    # Find products that are completed before t
    ytrain <- y[1:j]
    n_prior <- which((end < start[i]+days(7*(j-1))))
    if (length(n_prior) < 2) next # Move to the next one if there are less than two completed life cycles
    # Find the model corresponding to the set of completed life cycles
    n_model <-  which(lapply(priors_list_unique,function(x) length(x) == length(n_prior) && all(x == n_prior)) == TRUE)
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

# Print out average accuracy results
result_vec <- t(as.matrix(2*c(mean(rowMeans(apply(accuracy_result[,,1:8,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                  mean(rowMeans(apply(apply(accuracy_result[,,1:8,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                  mean(rowMeans(apply(accuracy_result[,,9:17,50],c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE),
                  mean(rowMeans(apply(apply(accuracy_result[,,9:17,1:99],c(1,2,3),mean),c(1,3),mean,na.rm=TRUE), na.rm=TRUE),na.rm=TRUE))))
colnames(result_vec) <- c('1-8 steps MAE','1-8 steps CRPS','9-17 steps MAE','9-17 steps CRPS')
result_vec


