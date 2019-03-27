library(smcfcs)
library(survival)
context("Weibull model testing")

test_that("Weibull imputation runs, binary covariate", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1*(runif(n)<0.5)
    z <- x+rnorm(n)
    t <- rsurvreg(n, mean=x+z, scale=2)
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x,z)

    imps <- smcfcs(simData, smtype="weibull", smformula="Surv(t, d)~x+z",
                   method=c("", "", "logreg", ""), m=1)
  }, NA)
})

test_that("Weibull imputation is consistent, binary covariate", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 100000
    x <- 1*(runif(n)<0.5)
    z <- x+rnorm(n)
    t <- rsurvreg(n, mean=x+z, scale=2)
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x,z)

    imps <- smcfcs(simData, smtype="weibull", smformula="Surv(t, d)~x+z",
                   method=c("", "", "logreg", ""), m=1)
    fitmod <- survreg(Surv(t,d)~x+z, data=imps$impDatasets[[1]])
    as.logical((abs(coef(fitmod)[2]-1)<0.05) & (abs(fitmod$scale-2)<0.05))
  }, TRUE)
})

test_that("Weibull imputation is consistent, cts covariate", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 100000
    x <- rnorm(n)
    z <- x+rnorm(n)
    t <- rsurvreg(n, mean=x+z, scale=2)
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.01)] <- NA

    simData <- data.frame(t,d,x,z)

    imps <- smcfcs(simData, smtype="weibull", smformula="Surv(t, d)~x+z",
                   method=c("", "", "norm", ""), m=1)
    fitmod <- survreg(Surv(t,d)~x+z, data=imps$impDatasets[[1]])
    as.logical((abs(coef(fitmod)[2]-1)<0.05) & (abs(fitmod$scale-2)<0.05))
  }, TRUE)
})
