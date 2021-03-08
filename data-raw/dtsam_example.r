#the following example is not run when the package is compiled on CRAN
#(to keep computation time down), but it can be run by package users
\dontrun{
  #discrete time survival analysis example
  M <- 5
  imps <- smcfcs.dtsam(ex_dtsam, "Surv(failtime,d)~x1+x2",
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
  library(mitools)
  summary(MIcombine(ests,vars))
}
