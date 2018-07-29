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

test_that("Checking error trap for method statement", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- y <- x1+x2+rnorm(n)
    x1[(runif(n)<0.5)] <- NA

    simData <- data.frame(x1,x2,y)

    imps <- smcfcs(simData, smtype="lm", smformula="y~x1+x2",
                   method=c("", "", ""))
  })
})

test_that("Checking error trap for method statement", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- y <- x1+x2+rnorm(n)
    x1[(runif(n)<0.5)] <- NA

    simData <- data.frame(x1,x2,y)

    imps <- smcfcs(simData, smtype="lm", smformula="y~x1+x2",
                   method=c("", "norm", ""))
  })
})

test_that("Checking error trap for method statement", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- y <- x1+x2+rnorm(n)
    x1[(runif(n)<0.5)] <- NA

    simData <- data.frame(x1,x2,y)

    imps <- smcfcs(simData, smtype="lm", smformula="y~x1+x2",
                   method=c("norm", "", ""))
  })
})
