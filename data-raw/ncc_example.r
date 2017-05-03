#the following example is not run when the package is compiled on CRAN
#(to keep computation time down), but it can be run by package users
\dontrun{
  predictorMatrix <- matrix(0,nrow=dim(ex_ncc)[2],ncol=dim(ex_ncc)[2])
  predictorMatrix[which(colnames(ex_ncc)=="x"),c(which(colnames(ex_ncc)=="z"))] <- 1

  imps <- smcfcs.nestedcc(originaldata=ex_ncc,set="setno",nrisk="numrisk",event="d",
                          smformula="Surv(t,case)~x+z+strata(setno)",
                          method=c("", "", "logreg", "", "", "", "", ""),
                          predictorMatrix=predictorMatrix)
  library(mitools)
  impobj <- imputationList(imps$impDatasets)
  models <- with(impobj, clogit(case~x+z+strata(setno)))
  summary(MIcombine(models))
}
