#the following example is not run when the package is compiled on CRAN
#(to keep computation time down), but it can be run by package users
\dontrun{
  #as per the documentation for ex_cc, the sampling fraction is 10%
  imps <- smcfcs.casecohort(ex_cc, smformula="Surv(entertime, t, d)~x+z", sampfrac=0.1,
                            in.subco="in.subco", method=c("", "", "norm", "", "", "", ""))
  library(mitools)
  impobj <- imputationList(imps$impDatasets)
  models <- with(impobj, coxph(Surv(entertime,t,d)~x+z+cluster(id)))
  summary(MIcombine(models))
}
