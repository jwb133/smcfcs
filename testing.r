imps <- smcfcs.lm(ex_linquad, smformula="y~z+x+xsq",method=c("","","norm","x^2"),m=10)

library(mitools)
impobj <- imputationList(imps)
models <- with(impobj, lm(y~z+x+xsq))
summary(MIcombine(models))


imps <- smcfcs.lm(ex_lininter, smformula="y~x1+x2+x1x2",method=c("","norm","logreg","x1*x2"),m=10)

library(mitools)
impobj <- imputationList(imps)
models <- with(impobj, lm(y~x1+x2+x1x2))
summary(MIcombine(models))
