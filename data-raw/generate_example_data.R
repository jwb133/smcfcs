# linear substantive model with quadratic covariate effect
set.seed(1234)
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
# v is an auxiliary variable
v <- x + rnorm(n)
y <- 1 + z + x + x^2 + rnorm(n)

# make some x values missing
xobsxb <- (y - mean(y)) / sd(y)
xobspr <- exp(xobsxb) / (1 + exp(xobsxb))
x[runif(n) > xobspr] <- NA

ex_linquad <- data.frame(y, z, x, v)

usethis::use_data(ex_linquad, overwrite = TRUE)

# linear substantive model with interaction
n <- 1000
x1 <- rnorm(n)
x2 <- 1 * (runif(n) < 0.5)
y <- 1 + x1 + x2 + x1 * x2 + rnorm(n)

# make some x1 and x2 values missing
x1[runif(n) > 0.5] <- NA
x2[runif(n) > 0.5] <- NA

ex_lininter <- data.frame(y, x1, x2)
ex_lininter$x2 <- factor(ex_lininter$x2)

usethis::use_data(ex_lininter, overwrite = TRUE)

# logistic substantive model with quadratic covariate effect
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
# v is an auxiliary variable
v <- x + rnorm(n)
xb <- z + x + x^2
pr <- exp(xb) / (1 + exp(xb))
y <- 1 * (runif(n) < pr)

# make some x values missing
xobsxb <- (z - mean(z)) / sd(z)
xobspr <- exp(xobsxb) / (1 + exp(xobsxb))
x[runif(n) > xobspr] <- NA

ex_logisticquad <- data.frame(y, z, x, v)

usethis::use_data(ex_logisticquad, overwrite = TRUE)

# Cox substantive model with quadratic covariate effect
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
# v is an auxiliary variable
v <- x + rnorm(n)
xb <- z + x + x^2
c <- -log(runif(n)) / 1
t <- -log(runif(n)) / exp(xb)
d <- 1 * (t < c)
t[d == 0] <- c[d == 0]

# make some x values missing
xobsxb <- (t - mean(t)) / sd(t)
xobspr <- exp(xobsxb) / (1 + exp(xobsxb))
x[runif(n) > xobspr] <- NA

ex_coxquad <- data.frame(t, d, z, x, v)

usethis::use_data(ex_coxquad, overwrite = TRUE)

# competing risks data
n <- 1000
x1 <- rbinom(n, 1, 0.5)
x2 <- x1 + rnorm(n)
xb1 <- x1 + x2
xb2 <- -x1
t1 <- -log(runif(n)) / (0.002 * exp(xb1))
t2 <- -log(runif(n)) / (0.002 * exp(xb2))
c <- -log(runif(n)) / 0.002
d <- rep(0, n)
d[(t1 < c) & (t1 < t2)] <- 1
d[(t2 < c) & (t2 < t1)] <- 2
t <- c
t[d == 1] <- t1[d == 1]
t[d == 2] <- t2[d == 2]
x1[runif(n) > 0.5] <- NA
x2[runif(n) > 0.5] <- NA

ex_compet <- data.frame(t, d, x1, x2)

usethis::use_data(ex_compet, overwrite = TRUE)

# Poisson regression
n <- 1000
x <- rnorm(n)
z <- 1 * (runif(n) < 0.5)
y <- rpois(n, lambda = exp(x - z))

# make some x1 and x2 values missing
x[runif(n) > (exp(y) / (1 + exp(y)))] <- NA

ex_poisson <- data.frame(y, x, z)

usethis::use_data(ex_poisson, overwrite = TRUE)

# case cohort
n <- 10000
z <- rnorm(n)
x <- z + rnorm(n)
t <- -log(runif(n)) / (0.01 * exp(x + z))
d <- 1 * (t < 1)
t[d == 0] <- 1
x[(runif(n) < 0.5)] <- NA

fullcohortdata <- data.frame(t, d, x, z)
fullcohortdata$in.subco <- 0
# we sample a 10% subcohort
fullcohortdata$in.subco[sample(n, size = n * 0.1)] <- 1
fullcohortdata$id <- 1:n

ex_cc <- fullcohortdata[(fullcohortdata$in.subco == 1) | (fullcohortdata$d == 1), ]
ex_cc$entertime <- 0
ex_cc$entertime[ex_cc$in.subco == 0] <- ex_cc$t[ex_cc$in.subco == 0] - 0.000001

usethis::use_data(ex_cc, overwrite = TRUE)

# nested case control
n <- 10000
z <- rnorm(n)
x <- rbinom(n, 1, exp(z) / (1 + exp(z)))
t <- -log(runif(n)) / (0.01 * exp(x + z))
d <- 1 * (t < 1)
t[d == 0] <- 1
x[(runif(n) < 0.5)] <- NA

fullcohortdata <- data.frame(t, d, x, z)
fullcohortdata$id <- 1:n

# Compute number at risk at each event time using the full cohort data
nrisk.fit <- survival::survfit(Surv(t, d) ~ 1, data = fullcohortdata)
ord.t.d1 <- order(fullcohortdata$t[fullcohortdata$d == 1])

m <- 1 # 1 control per case
ex_ncc <- NULL

no.sample <- 0
for (i in which(fullcohortdata$d == 1))
{
  # select control(s) for nested case-control
  possible.controls <- which(fullcohortdata$t >= fullcohortdata$t[i])
  # remove the case from this vector
  possible.controls <- possible.controls[which(possible.controls != i)]

  if (length(possible.controls) >= m) {
    controls <- sample(possible.controls, m)
    numAtRisk <- 1 + length(possible.controls)
    ex_ncc <- rbind(ex_ncc, cbind(fullcohortdata[i, ], numrisk = numAtRisk))
    ex_ncc <- rbind(ex_ncc, cbind(fullcohortdata[controls, ], numrisk = numAtRisk))
    no.sample <- no.sample + 1
  }
}

ex_ncc$setno <- rep(1:no.sample, each = m + 1)
ex_ncc$case <- rep(c(1, rep(0, m)), no.sample)

usethis::use_data(ex_ncc, overwrite = TRUE)

# discrete time survival analysis
n <- 1000
x1 <- 1 * (runif(n) < 0.5)
x2 <- x1 + rnorm(n)
# define number of time points
T <- 10
# define vector of intercepts
alpha <- seq(-1, -0.1, 0.1)
beta <- c(1, -1)
yMat <- array(0, dim = c(n, T))
for (i in 1:T) {
  yMat[, i] <- 1 * (runif(n) < expit(alpha[i] + beta[1] * x1 + beta[2] * x2))
}
failtime <- apply(yMat, 1, function(x) which(x == 1)[1])
# event indicator
d <- rep(1, n)
d[is.na(failtime)] <- 0
failtime[is.na(failtime)] <- 10
mean(d)
# add in some additional completely random censoring
cTime <- ceiling(10 * runif(n))
t <- failtime
t[cTime < failtime] <- cTime[cTime < failtime]
d[cTime < failtime] <- 0
mean(d)

# make some data missing in x1
x1[runif(n) < expit((x2 - mean(x2) / sd(x2)))] <- NA

ex_dtsam <- data.frame(x1 = x1, x2 = x2, failtime = t, d = d)
usethis::use_data(ex_dtsam, overwrite = TRUE)


# Fine-Gray example data --------------------------------------------------


n <- 1000

# Covariates the same as ex_compet
x1 <- rbinom(n = n, size = 1, prob = 0.5)
x2 <- rnorm(n = n, mean = x1)

# p = proportion cause 1 failures when x1 and x2 are zero
p <- 0.15
xb1 <- 0.75 * x1 - 0.5 * x2 # Fine-Gray linear predictor
xb2 <- -x1 + 0.5 * x2 # Linear predictor for hazard conditional on cause 2

# Generate competing event indicator
d_tilde <- 1 + rbinom(n = n, size = 1, prob = (1 - p)^exp(xb1))
cause1_ind <- d_tilde == 1

# Generate event times conditional on the indicator
U <- runif(n = n)
t_tilde <- numeric(length = n)
num <- 1 - (1 - U[cause1_ind] * (1 - (1 - p)^xb1[cause1_ind]))^(1 / xb1[cause1_ind])
t_tilde[cause1_ind] <- -log(1 - num / p)
t_tilde[!cause1_ind] <- -log(U[!cause1_ind]) / exp(xb2[!cause1_ind])
cens <- -log(runif(n = n)) / 0.1
times <- pmin(cens, t_tilde)
d <- ifelse(cens < t_tilde, 0, d_tilde)

# Same MCAR missingness as ex_compet
x1[runif(n) > 0.5] <- NA
x2[runif(n) > 0.5] <- NA

ex_finegray <- data.frame(times, d, x1, x2)
usethis::use_data(ex_finegray, overwrite = TRUE)


# Covariate measurement error ---------------------------------------------


# n <- 10000
# x <- rnorm(n)
# xb <- x
# pr <- exp(xb)/(1+exp(xb))
# y <- 1*(runif(n)<pr)
#
# w1 <- x+rnorm(n)
# w2 <- x+rnorm(n)
# w2[runif(1000)<0.9] <- NA
#
# ex_coverr <- data.frame(y,x=NA,w1,w2)
#
# usethis::use_data(ex_coverr, overwrite=TRUE)

# Flexsurv example data --------------------------------------------------

set.seed(123)
n <- 1000
myshape <- 1.5
myscale <- 3
z <- rnorm(n)
x <- 1*(runif(n)<expit(z))
t <- rweibull(n=n,shape=myshape,scale=exp(x+z)^(-1/myshape))
c <- rexp(n)
d <- 1*(t<c)
t[d==0] <- c[d==0]
x[runif(n)<0.5] <- NA

ex_flexsurv <- data.frame(t=t,d=d,x=x,z=z)
usethis::use_data(ex_flexsurv, overwrite = TRUE)


# Restrictions example data --------------------------------------------------

set.seed(723423)
n <- 100
z <- rnorm(n)

p1 <- 1/(1+exp(0+z)+exp(1+z))
p2 <- exp(0+z)/(1+exp(0+z)+exp(1+z))
p3 <- 1-p1-p2
u <- runif(n)
x <- rep("a",n)
x[(u>p1) & (u<(p1+p2))] <- "b"
x[(u>(p1+p2))] <- "c"

y <- z+(x=="b")+2*(x=="c")+rnorm(n)

# generate coarsening variable, where 0 indicates no coarsening
# 1 indicates missing
# 2 indicates coarsened variable=a if x=a or b/c if x=b or x=c
# 3 indicates coarsened variable=b if x=b or a/c if x=a or x=c
c <- sample(c(0, 1, 2, 3), size = n, replace = TRUE, prob = rep(0.25,4))
xobs <- x
xobs[c==1] <- "NA"
xobs[(c==2) & ((x=="b") | (x=="c"))] <- "b/c"
xobs[(c==3) & ((x=="a") | (x=="c"))] <- "a/c"
table(xobs, useNA = "always")

# set x to NA when not exactly observed
x[c==1] <- NA
x[(c==2) & ((x=="b") | (x=="c"))] <- NA
x[(c==3) & ((x=="a") | (x=="c"))] <- NA

ex_coarsening <- data.frame(x=factor(x),xobs=xobs,z=z,y=y)
usethis::use_data(ex_coarsening, overwrite = TRUE)

