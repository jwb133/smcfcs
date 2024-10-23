library(smcfcs)
library(survival)
context("Flexible parametric proportional hazards model testing")

test_that("Flexsurv imputation of missing normal covariate is approximately unbiased", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(1234)
      n <- 1000
      z <- rnorm(n)
      x <- z + rnorm(n)
      t <- -log(runif(n)) / (1 * exp(x + z))
      d <- 1 * (t < 10)
      t[d == 0] <- 10
      x[(runif(n) < 0.5)] <- NA

      simData <- data.frame(t, d, x, z)

      imps <- smcfcs.flexsurv(simData,
        smformula = "Surv(t, d)~x+z",
        method = c("", "", "norm", ""),
        k=2
      )
      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, flexsurvspline(Surv(t, d) ~ x + z, k=2))
      abs(summary(MIcombine(models))[5, 1] - 1) < 0.1
    },
    TRUE
  )

})

test_that("Flexsurv imputation of missing binary covariate is approximately unbiased", {
  skip_on_cran()
  expect_equal(
    {
      expit <- function(x) {exp(x)/(1+exp(x))}
      set.seed(1234)
      n <- 1000
      z <- rnorm(n)
      x <- 1*(runif(n)<expit(z))
      t <- -log(runif(n)) / (1 * exp(x + z))
      d <- 1 * (t < 10)
      t[d == 0] <- 10
      x[(runif(n) < 0.5)] <- NA

      simData <- data.frame(t, d, x, z)

      imps <- smcfcs.flexsurv(simData,
                              smformula = "Surv(t, d)~x+z",
                              method = c("", "", "logreg", ""),
                              k=2
      )
      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, flexsurvspline(Surv(t, d) ~ x + z, k=2))
      abs(summary(MIcombine(models))[5, 1] - 1) < 0.1
    },
    TRUE
  )

})



