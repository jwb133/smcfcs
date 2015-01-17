#linear substantive models
imps <- smcfcs.lm(ex_linquad, smformula="y~z+x+xsq",method=c("","","norm","x^2",""),m=10)

library(mitools)
impobj <- imputationList(imps)
models <- with(impobj, lm(y~z+x+xsq))
summary(MIcombine(models))

#include auxiliary variable assuming it is conditionally independent of Y (which it is here)
predMatrix <- array(0, dim=c(ncol(ex_linquad),ncol(ex_linquad)))
predMatrix[3,] <- c(0,1,0,0,1)
imps <- smcfcs.lm(ex_linquad, smformula="y~z+x+xsq",method=c("","","norm","x^2",""),predictorMatrix=predMatrix,m=10)

#interaction model
imps <- smcfcs.lm(ex_lininter, smformula="y~x1+x2+x1x2",method=c("","norm","logreg","x1*x2"),m=10)

library(mitools)
impobj <- imputationList(imps)
models <- with(impobj, lm(y~x1+x2+x1x2))
summary(MIcombine(models))
