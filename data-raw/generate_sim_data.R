#linear substantive model with quadratic covariate effect
set.seed(1234)
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
y <- 1+z+x+x^2+rnorm(n)

#make some x values missing
xobsxb <- (y-mean(y))/sd(y)
xobspr <- exp(xobsxb)/(1+exp(xobsxb))
x[runif(n)>xobspr] <- NA

ex_linquad <- data.frame(y,z,x,xsq=x^2)

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
