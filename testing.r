#linear substantive models
imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+xsq",method=c("","","norm","x^2",""))

library(mitools)
impobj <- imputationList(imps)
models <- with(impobj, lm(y~z+x+xsq))
summary(MIcombine(models))

#include auxiliary variable assuming it is conditionally independent of Y (which it is here)
predMatrix <- array(0, dim=c(ncol(ex_linquad),ncol(ex_linquad)))
predMatrix[3,] <- c(0,1,0,0,1)
imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+xsq",method=c("","","norm","x^2",""),predictorMatrix=predMatrix)

#interaction model
imps <- smcfcs(ex_lininter, smtype="lm", smformula="y~x1+x2+x1x2",method=c("","norm","logreg","x1*x2"))

library(mitools)
impobj <- imputationList(imps)
models <- with(impobj, lm(y~x1+x2+x1x2))
summary(MIcombine(models))

#logistic regression substantive model
imps <- smcfcs(ex_logisticquad, smtype="logistic", smformula="y~z+x+xsq",method=c("","","norm","x^2",""))

library(mitools)
impobj <- imputationList(imps)
models <- with(impobj, glm(y~z+x+xsq, family=binomial))
summary(MIcombine(models))

#Cox regression substantive model
imps <- smcfcs(ex_coxquad, smtype="coxph", smformula="Surv(t,delta)~z+x+xsq",method=c("","","","norm","x^2",""))

library(mitools)
impobj <- imputationList(imps)
models <- with(impobj, coxph(Surv(t,delta)~z+x+xsq))
summary(MIcombine(models))
