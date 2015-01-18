#linear substantive model with quadratic covariate effect
set.seed(1234)
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
#v is an auxiliary variable
v <- x+rnorm(n)
y <- 1+z+x+x^2+rnorm(n)

#make some x values missing
xobsxb <- (y-mean(y))/sd(y)
xobspr <- exp(xobsxb)/(1+exp(xobsxb))
x[runif(n)>xobspr] <- NA

ex_linquad <- data.frame(y,z,x,xsq=x^2,v)

devtools::use_data(ex_linquad, overwrite=TRUE)

#linear substantive model with interaction
set.seed(1234)
n <- 1000
x1 <- rnorm(n)
x2 <- 1*(runif(n)<0.5)
y <- 1+x1+x2+x1*x2+rnorm(n)

#make some x1 and x2 values missing
x1[runif(n)>0.5] <- NA
x2[runif(n)>0.5] <- NA

ex_lininter <- data.frame(y,x1,x2,x1x2=x1*x2)

devtools::use_data(ex_lininter, overwrite=TRUE)

#logistic substantive model with quadratic covariate effect
set.seed(1234)
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
#v is an auxiliary variable
v <- x+rnorm(n)
xb <- z+x+x^2
pr <- exp(xb)/(1+exp(xb))
y <- 1*(runif(n)<pr)

#make some x values missing
xobsxb <- (z-mean(z))/sd(z)
xobspr <- exp(xobsxb)/(1+exp(xobsxb))
x[runif(n)>xobspr] <- NA

ex_logisticquad <- data.frame(y,z,x,xsq=x^2,v)

devtools::use_data(ex_logisticquad, overwrite=TRUE)

#Cox substantive model with quadratic covariate effect
#set.seed(1234)
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
#v is an auxiliary variable
v <- x+rnorm(n)
xb <- z+x+x^2
t <- -log(runif(n))/exp(xb)
delta <- 1*(t<20)
t[delta==0] <- 20

#make some x values missing
xobsxb <- (t-mean(t))/sd(t)
xobspr <- exp(xobsxb)/(1+exp(xobsxb))
x[runif(n)>xobspr] <- NA

ex_coxquad <- data.frame(t,delta,z,x,xsq=x^2,v)

devtools::use_data(ex_coxquad, overwrite=TRUE)

#competing risks data
n <- 1000
x1 <- rbinom(n,1,0.5)
x2 <- x1+rnorm(n)
xb1 <- x1+x2
xb2 <- -x1
t1 <- -log(runif(n))/(0.002*exp(xb1))
t2 <- -log(runif(n))/(0.002*exp(xb2))
c <- -log(runif(n))/0.002
d <- rep(0,n)
d[(t1<c) & (t1<t2)] <- 1
d[(t2<c) & (t2<t1)] <- 2
t <- c
t[d==1] <- t1[d==1]
t[d==2] <- t2[d==2]
x1[runif(n)>0.5] <- NA
x2[runif(n)>0.5] <- NA

ex_compet <- data.frame(t,d,x1,x2)

devtools::use_data(ex_compet, overwrite=TRUE)
