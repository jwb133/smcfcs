#set random number seed to make results reproducible
set.seed(123)

#simulate a dataset
n <- 200
x <- rnorm(n)
psxb <- x
pspr <- exp(psxb)/(1+exp(psxb))
z <- 1*(runif(n)<pspr)
y <- x + z + rnorm(n)
data <- data.frame(x,z,y)

#make some x values missing completely at random
data$x[runif(n)<0.5] <- NA

#generate multiple imputations using mips
imps <- mips(data, omtype="lm", omformula="y~x+z",
             method=c("norm","",""),psformula="z~x",m=5)
#now do some propensity score analysis using the imputed datasets
