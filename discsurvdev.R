n <- 100000
library(MASS)

#source Caro's file defining genDiscSurv - this is a slightly modified version
#where I got rid of the bit which modified the first 30 persons to avoid perfect
#prediction (I did this just in case it was the source of the slight full data biases
#I think I was seeing)
source("genDiscsurv.R")

simData <- data.frame(genDiscSurv(n=n))
head(simData)
summary(simData)

longFullData <- survSplit(Surv(y,failed)~X1+X2+X3+X4+X5, data=simData, cut=unique(simData$y[simData$failed==1]))
mod <- glm(failed~-1+factor(tstart)+X1+X2+X3+X4+X5, family="binomial", data=longFullData)
summary(mod)

#make some data missing completely at random
simData$X1[runif(n)<0.5] <- NA
simData$X2[runif(n)<0.5] <- NA
simData$X3[runif(n)<0.5] <- NA
simData$X4[runif(n)<0.5] <- NA

library(survival)
imps <- smcfcs(simData, "dtsam", "Surv(y,failed)~X1+X2+X3+X4+X5",
               method=c(rep("norm",4),rep("",4)),m=1)

#fit outcome model to first imputed dataset. Need to code up fitting using all of them
#and combining using Rubin's rules
#reshape to long form
longData <- survSplit(Surv(y,failed)~X1+X2+X3+X4+X5, data=imps$impDatasets[[1]], cut=unique(simData$y[simData$failed==1]))
mod <- glm(failed~-1+factor(tstart)+X1+X2+X3+X4+X5, family="binomial", data=longData)
summary(mod)
