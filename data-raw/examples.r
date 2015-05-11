#load mitools for fitting models to imputed datasets
library(mitools)

#set random number seed to make results reproducible
set.seed(123)

#linear substantive model with quadratic covariate effect
imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+xsq",
               method=c("","","norm","x^2",""))
impobj <- imputationList(imps$impDatasets)
models <- with(impobj, lm(y~z+x+xsq))
summary(MIcombine(models))

#examining convergence, using 100 iterations, setting m=1
imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+xsq",
               method=c("","","norm","x^2",""),m=1,numit=100)
#convergence plot from first imputation for third coefficient of substantive model
plot(imps$smCoefIter[1,3,])

#include auxiliary variable assuming it is conditionally independent of Y (which it is here)
predMatrix <- array(0, dim=c(ncol(ex_linquad),ncol(ex_linquad)))
predMatrix[3,] <- c(0,1,0,0,1)
imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+xsq",
               method=c("","","norm","x^2",""),predictorMatrix=predMatrix)

#interaction model
imps <- smcfcs(ex_lininter, smtype="lm", smformula="y~x1+x2+x1x2",
               method=c("","norm","logreg","x1*x2"))
impobj <- imputationList(imps$impDatasets)
models <- with(impobj, lm(y~x1+x2+x1x2))
summary(MIcombine(models))

#logistic regression substantive model
imps <- smcfcs(ex_logisticquad, smtype="logistic", smformula="y~z+x+xsq",
               method=c("","","norm","x^2",""))
impobj <- imputationList(imps$impDatasets)
models <- with(impobj, glm(y~z+x+xsq, family=binomial))
summary(MIcombine(models))

#Cox regression substantive model
library(survival)
imps <- smcfcs(ex_coxquad, smtype="coxph", smformula="Surv(t,d)~z+x+xsq",
               method=c("","","","norm","x^2",""))
impobj <- imputationList(imps$impDatasets)
models <- with(impobj, coxph(Surv(t,d)~z+x+xsq))
summary(MIcombine(models))

#competing risks substantive model
imps <- smcfcs(ex_compet, smtype="compet", smformula=c("Surv(t,d==1)~x1+x2", "Surv(t,d==2)~x1+x2"),
               method=c("","","logreg","norm"))
impobj <- imputationList(imps$impDatasets)
models <- with(impobj, coxph(Surv(t,d==1)~x1+x2))
summary(MIcombine(models))
models <- with(impobj, coxph(Surv(t,d==2)~x1+x2))
summary(MIcombine(models))
