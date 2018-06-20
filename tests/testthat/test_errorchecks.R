library(smcfcs)
library(survival)
context("Error trap testing")

test_that("Checking outcome model check for logistic models", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- as.factor(1*(runif(n)<0.5))
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(x,y)

    imps <- smcfcs(simData, smtype="logistic", smformula="y~x",
                   method=c("norm", ""))
  })
})

test_that("Checking outcome model check for logistic models", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- 1+1*(runif(n)<0.5)
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(x,y)

    imps <- smcfcs(simData, smtype="logistic", smformula="y~x",
                   method=c("norm", ""))
  })
})

test_that("Checking outcome model check for logistic models", {
  expect_output({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- 1*(runif(n)<0.5)
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(x,y)

    imps <- smcfcs(simData, smtype="logistic", smformula="y~x",
                   method=c("norm", ""))
  })
})
