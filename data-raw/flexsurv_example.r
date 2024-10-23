#the following example is not run when the package is compiled on CRAN
#(to keep computation time down), but it can be run by package users
\dontrun{

  set.seed(63213)
  imps <- smcfcs.flexsurv(ex_flexsurv,
                          k=2,
                          smformula="Surv(t,d)~x+z",
                          method=c("","","logreg",""))
  library(mitools)
  impobj <- imputationList(imps$impDatasets)
  models <- with(impobj, flexsurvspline(Surv(t,d)~x+z, k=2))
  summary(MIcombine(models))
}
