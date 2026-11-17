library(smcfcs)
library(mitools)
library(survival)
context("Poisson model testing")

test_that("smcfcs runs with missing Poisson outcomes", {
  expect_error(
    {
      set.seed(1234)
      n <- 100
      x <- 1 * (runif(n) < 0.5)
      z <- x + rnorm(n)
      y <- rpois(n=n, lambda=exp(x+z))
      y[1:50] <- NA
      x[1:10] <- NA

      simData <- data.frame(y, x, z)

      imps <- smcfcs(simData,
        smtype = "poisson", smformula = "y~x+z",
        method = c("", "logreg", ""), m = 1
      )
    },
    NA
  )
})

test_that("smcfcs is unbiased with missing Poisson outcomes", {
  skip_on_cran()
  expect_error(
    {
      set.seed(1234)
      n <- 10000
      x <- 1 * (runif(n) < 0.5)
      z <- x + rnorm(n)
      y <- rpois(n=n, lambda=exp(x+z))
      y[1:5000] <- NA
      x[1:10] <- NA

      simData <- data.frame(y, x, z)

      imps <- smcfcs(simData,
                     smtype = "poisson", smformula = "y~x+z",
                     method = c("", "logreg", ""), m = 10
      )
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, glm(y~x+z,family=poisson))
      res <- summary(MIcombine(models))
      as.logical((res$`(lower`[2] < 1) & (res$`upper)`[2] < 1))
    },
    NA
  )
})
