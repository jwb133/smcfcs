#set random number seed to make results reproducible
set.seed(123)

#linear substantive model with quadratic covariate effect
imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+I(x^2)",
               method=c("","","norm",""))

#if mitools is installed, fit substantive model to imputed datasets
#and combine results using Rubin's rules
if (requireNamespace("mitools", quietly = TRUE)) {
  library(mitools)
  impobj <- imputationList(imps$impDatasets)
  models <- with(impobj, lm(y~z+x+I(x^2)))
  summary(MIcombine(models))
}

#the following examples are not run when the package is compiled on CRAN
#(to keep computation time down), but they can be run by package users
\dontrun{
  #examining convergence, using 100 iterations, setting m=1
  imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+I(x^2)",
                 method=c("","","norm",""),m=1,numit=100)
  #convergence plot from first imputation for third coefficient of substantive model
  plot(imps$smCoefIter[1,3,])

  #include auxiliary variable assuming it is conditionally independent of Y (which it is here)
  predMatrix <- array(0, dim=c(ncol(ex_linquad),ncol(ex_linquad)))
  predMatrix[3,] <- c(0,1,0,1)
  imps <- smcfcs(ex_linquad, smtype="lm", smformula="y~z+x+I(x^2)",
                 method=c("","","norm",""),predictorMatrix=predMatrix)

  #impute missing x1 and x2, where they interact in substantive model
  imps <- smcfcs(ex_lininter, smtype="lm", smformula="y~x1+x2+x1*x2",
                 method=c("","norm","logreg"))

  #logistic regression substantive model, with quadratic covariate effects
  imps <- smcfcs(ex_logisticquad, smtype="logistic", smformula="y~z+x+xsq",
                 method=c("","","norm","x^2",""))

  #Poisson regression substantive model
  imps <- smcfcs(ex_poisson, smtype="poisson", smformula="y~x+z",
                 method=c("","norm",""))
  if (requireNamespace("mitools", quietly = TRUE)) {
    library(mitools)
    impobj <- imputationList(imps$impDatasets)
    models <- with(impobj, glm(y~x+z,family=poisson))
    summary(MIcombine(models))
  }

  #Cox regression substantive model, with only main covariate effects
  if (requireNamespace("survival", quietly = TRUE)) {
    imps <- smcfcs(ex_coxquad, smtype="coxph", smformula="Surv(t,d)~z+x+xsq",
                   method=c("","","","norm","x^2",""))

    #competing risks substantive model, with only main covariate effects
    imps <- smcfcs(ex_compet, smtype="compet",
                   smformula=c("Surv(t,d==1)~x1+x2", "Surv(t,d==2)~x1+x2"),
                   method=c("","","logreg","norm"))
  }

  #if mitools is installed, fit model for first competing risk
  if (requireNamespace("mitools", quietly = TRUE)) {
    library(mitools)
    impobj <- imputationList(imps$impDatasets)
    models <- with(impobj, coxph(Surv(t,d==1)~x1+x2))
    summary(MIcombine(models))
  }

  #discrete time survival analysis example
  M <- 5
  imps <- smcfcs(ex_dtsam, "dtsam", "Surv(failtime,d)~x1+x2",
                 method=c("logreg","", "", ""),m=M)
  #fit dtsam model to each dataset manually, since we need
  #to expand to person-period data form first
  ests <- vector(mode = "list", length = M)
  vars <- vector(mode = "list", length = M)
  for (i in 1:M) {
    longData <- survSplit(Surv(failtime,d)~x1+x2, data=imps$impDatasets[[i]],
                          cut=unique(ex_dtsam$failtime[ex_dtsam$d==1]))
    mod <- glm(d~-1+factor(tstart)+x1+x2, family="binomial", data=longData)
    ests[[i]] <- coef(mod)
    vars[[i]] <- diag(vcov(mod))
  }
  summary(MIcombine(ests,vars))

}
