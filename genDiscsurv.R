
## Filename: genDiscsurv.R
## Description: Create independent variables and
## create discrete survival data for simulation
### This function WITHOUT ADDED FRAILTY

# n Number of observations for the sample
# r Correlation between independent variables
# P_max Number of maximal possible time points after which
# all observeations are censored
# censoring Proportion of observations that are randomly censored




## R version 3.5.3 (2019-03-11)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 7 x64 (build 7601) Service Pack 1

## Matrix products: default

## Creation: 2019-03-12


genDiscSurv <- function(n = numberobs,
                           r = 0.1,
                           P_max =16,
                           censoring=0.1){

# Create indizes matrices
P_ind_single <- matrix(rep(rep(0,P_max),P_max),nrow = P_max)
diag(x = P_ind_single) <- 1
P_ind_n <- do.call("rbind", replicate(n, P_ind_single, simplify = FALSE))

# Coefficients for the creation of survival times
alpha <- seq(from=-5,to=-0.8,length.out = P_max)
alpha <- matrix(alpha, ncol = 1)
beta <- c(0.8, 2.2, -0.5, 0.3, -1.4)
beta <- matrix(beta, ncol = 1)

# Creation of independent variables
r <- 0.1
p <- length(beta)
R <- matrix(c(rep(r, p^2)), ncol = p)
diag(R) <- 1
R
mu <- rep(1, p)
SD <- rep(1, p)
S <- R * (SD %*% t(SD))
X <- mvrnorm(n, mu, S)
X <- cbind(X,id=1:n)
# Transform to long format
X_long <- do.call("rbind", replicate(P_max, X, simplify = FALSE))

X_long <- X_long[order(X_long[,6]),]
X_long <- X_long[,-(p+1)]

# Calculation of hazards for individuals in each period
hi <- 1/
         (
           1 +
          exp(
            -(
              P_ind_n%*%alpha + X_long %*% beta
            )
          )
         )

# Determine when event happens
unif_event <- runif(min = 0, max = 1,n = n*P_max)
failed <- as.data.frame(matrix(hi>unif_event, ncol =P_max, byrow = T))

# Take first event and save as time-to-event
colnames(failed) <- 1:P_max
TTE <- as.numeric(colnames(failed)
                  [max.col(failed, "first")*P_max^!rowSums(failed)])

# Censor about 10% of the survival times
unif_censor <- runif(min = 0, max = 1,n = n*P_max)
censor_ind <- as.data.frame(matrix(unif_censor>(1-censoring),
                                   ncol =P_max, byrow = T))
colnames(censor_ind) <- 1:P_max
censor_p <- as.numeric(colnames(censor_ind)
           [max.col(censor_ind, "first")*P_max^!rowSums(censor_ind)])

# Bind everything together and return
failed <- TTE<=censor_p
y <- ifelse(failed, TTE, censor_p)
failed[TTE==P_max] <- 0
y[TTE==P_max] <- P_max-1

  #Perfect Predictors verhindern
  # X[1:30, 1] <- 0
  # X[1:30, 2] <- 0
  # X[1:30, 3] <- 0
  # X[1:30, 4] <- 0
  # X[1:30, 5] <- 0
  # y[1:30] <- rep(1:15,2)
  # failed[1:30] <- c(rep(1,15), rep(0,15))

X <- cbind(X,y)
X <- cbind(X,failed)
X <- as.data.frame(X)
X_final <- as.matrix(data.frame(X1=X$V1,
                      X2=X$V2,
                      X3=X$V3,
                      X4=X$V4,
                      X5=X$V5,
                      y=X$y,
                      failed=X$failed,
                      id=X$id))

return(X_final)
}
