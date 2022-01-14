# ----------------------------
# Data preparation in Study 2
# ----------------------------

# See Section A.7 for more detail.

# Read in the raw data and remove dates and SKU numbers. 
data1 <- read.csv('dell_data_raw.csv')
timeSeriesDell <- data1[,-c(1:3)]

# Calculate the length of each series
length_rep <- rep(NA,170)
for (j in 1:170) {
  y1 <- as.numeric(timeSeriesDell[j, ])
  y <- c(na.omit(y1))
  length_rep[j] <- length(y)
}

# Create a matrix of life cycle length based on Table 7. In addition to columns 2-4 in the table, we also add the minimum and maximum PLC lengths to the matrix.  
length_matrix <- t(matrix(c(48,53,85,47,58,75,34,49,62,23,54,61,22,39,58),3,5))
min_val <- min(length_rep)
max_val <- max(length_rep)
length_matrix <- cbind(rep(min_val,5),length_matrix,rep(max_val,5))
mean_weeklyorder <- c(1021,456,267,166,66)

# Unnormalize the product life cycle by sampling from its distribution of cumulative sales given length. 
weekly_order <- rep(NA,170)
set.seed(201)
for (j in 1:170) {
  y1 <- as.numeric(timeSeriesDell[j, ])
  y <- c(na.omit(y1))
  n_length <- length(y)
  
# Calculate the conditional distribution of sales given a product's life cycle length, Pr(Q=q|N=n), based on Eq. (25).
  den <- rep(NA,5)
# For each of the five quintiles of life cycle sales, calculate Pr(N=n|Q=q). We assume that the life cycle length follows a uniform distribution between two consecutive quartiles. 
  for (i in 1:5) { 
    n <- which(length_matrix[i,] >= n_length)[1] 
    if (n == 1) n <- 2
    den[i] <- 1/(length_matrix[i,n] - length_matrix[i,n-1])
  }
# Calculate Pr(Q=q|N=n)
  prob <- den/sum(den) 
  
# Sample a quintile from the conditional distribution of sales. 
  sample_quint <- sample(1:5,1,prob=prob)

# We assume that sales follows a normal distribution with mean given in Table 7 and variance equals 10% of the mean value.
  weekly_order[j] <- rnorm(1,mean=mean_weeklyorder[sample_quint],sd=mean_weeklyorder[sample_quint]*0.1)*n_length
}

write.csv(weekly_order,'weekly_M.csv',row.names = FALSE)
