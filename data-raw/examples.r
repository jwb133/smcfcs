#load mitools for fitting models to imputed datasets
library(mitools)

#linear substantive model with quadratic covariate effect
imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+xsq",method=c("","","norm","x^2",""))
impobj <- imputationList(imps)
models <- with(impobj, lm(y~z+x+xsq))
summary(MIcombine(models))

#include auxiliary variable v assuming it is conditionally independent of y
predMatrix <- array(0, dim=c(ncol(ex_linquad),ncol(ex_linquad)))
predMatrix[3,] <- c(0,1,0,0,1)
imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+xsq",method=c("","","norm","x^2",""),predictorMatrix=predMatrix)

#interaction model
imps <- smcfcs(ex_lininter, smtype="lm", smformula="y~x1+x2+x1x2",method=c("","norm","logreg","x1*x2"))
impobj <- imputationList(imps)
models <- with(impobj, lm(y~x1+x2+x1x2))
summary(MIcombine(models))

#logistic regression substantive model
imps <- smcfcs(ex_logisticquad, smtype="logistic", smformula="y~z+x+xsq",method=c("","","norm","x^2",""))
impobj <- imputationList(imps)
models <- with(impobj, glm(y~z+x+xsq, family=binomial))
summary(MIcombine(models))

#Cox regression substantive model
imps <- smcfcs(ex_coxquad, smtype="coxph", smformula="Surv(t,delta)~z+x+xsq",method=c("","","","norm","x^2",""))
impobj <- imputationList(imps)
models <- with(impobj, coxph(Surv(t,delta)~z+x+xsq))
summary(MIcombine(models))

#competing risks substantive model
imps <- smcfcs(ex_compet, smtype="compet", smformula=c("Surv(t,d==1)~x1+x2", "Surv(t,d==2)~x1+x2"),method=c("","","logreg","norm"))
impobj <- imputationList(imps)
models <- with(impobj, coxph(Surv(t,d==1)~x1+x2))
summary(MIcombine(models))
models <- with(impobj, coxph(Surv(t,d==2)~x1+x2))
summary(MIcombine(models))
