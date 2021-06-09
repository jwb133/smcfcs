library(smcfcs)
library(survival)
context("Testing parallel computation")

test_that("Cox imputation runs in parallel", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 10000
    z <- rnorm(n)
    x <- z+rnorm(n)
    t <- -log(runif(n))/(1*exp(x+z))
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x,z)

    imps <- smcfcs.parallel(originaldata=simData, smtype="coxph", smformula="Surv(t, d)~x+z",
                   method=c("", "", "norm", ""), m=5, n_cores=5)
    library(mitools)
    impobj <- imputationList(imps$impDatasets)
    models <- with(impobj, coxph(Surv(t,d)~x+z))
    abs(summary(MIcombine(models))[1,1]-1)<0.1
  }, TRUE)
})

test_that("Cox imputation works with only one covariate", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 10000
    x <- rnorm(n)
    t <- -log(runif(n))/(0.01*exp(x))
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x)

    imps <- smcfcs(simData, smtype="coxph", smformula="Surv(t, d)~x",
                   method=c("", "", "norm"))
    library(mitools)
    impobj <- imputationList(imps$impDatasets)
    models <- with(impobj, coxph(Surv(t,d)~x))
    abs(summary(MIcombine(models))[1,1]-1)<0.1
  }, TRUE)
})
