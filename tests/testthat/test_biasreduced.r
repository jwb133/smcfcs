library(smcfcs)
context("Bias reduced model testing")

test_that("smcfcs gives error with perfect prediction and logreg method", {
  #note that this is not guaranteed (for an arbitrary seed) to give an error
  #because smcfcs initialises by randomly selecting some observed values, and
  #these are used in the covariate model fits (unlike in mice)
  suppressWarnings(
  expect_error({
    n <- 10
    z <- c(0,0,0,0,0,1,1,1,1,1)
    x <- z
    y <- x+rnorm(n)
    x[c(5,10)] <- NA

    simData <- data.frame(x,z,y)
    #impute x using z as an auxiliary variable
    predMat <- array(0, dim=c(3,3))
    predMat[1,2] <- 1
    rm(x,z,y)
    set.seed(62377)
    imps <- smcfcs(simData, smtype="lm", smformula="y~x",
                   method=c("logreg", "", ""), predictorMatrix=predMat, m=1)
  }))
})

test_that("smcfcs gives no warnings or errors with perfect prediction and brlogreg method", {
  expect_error({
    n <- 10
    z <- c(0,0,0,0,0,1,1,1,1,1)
    x <- z
    y <- x+rnorm(n)
    x[c(5,10)] <- NA

    simData <- data.frame(x,z,y)
    rm(x,z,y)
    #impute x using z as an auxiliary variable
    predMat <- array(0, dim=c(3,3))
    predMat[1,2] <- 1
    set.seed(62377)
    imps <- smcfcs(simData, smtype="lm", smformula="y~x",
                   method=c("brlogreg", "", ""), predictorMatrix=predMat, m=1)
  }, NA)
})

test_that("brlogreg gives results with minimal bias under perfect prediction between predictors", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 1000
    nSim <- 100
    expit <- function(x) {
      exp(x)/(1+exp(x))
    }

    simEsts <- array(0, dim=c(nSim,2))

    for (i in 1:nSim) {
      z <- 1*(runif(n)<0.5)
      #true log odds parameters of x|z are finite but very large
      x <- 1*(runif(n)<expit(-10+20*z))
      x[runif(n)<0.25] <- NA
      y <- x+rnorm(n)

      simData <- data.frame(x,z,y)
      rm(x,z,y)
      #impute x using z as an auxiliary variable
      predMat <- array(0, dim=c(3,3))
      predMat[1,2] <- 1

      imps <- smcfcs(simData, smtype="lm", smformula="y~x",
                     method=c("brlogreg", "", ""), predictorMatrix=predMat)
      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      smcfcsRes <- MIcombine(with(impobj, lm(y~x)))
      simEsts[i,] <- smcfcsRes$coefficients
    }
    ttestRes <- t.test(simEsts[,2], mu=1)
    ttestRes$p.value>0.05
  }, TRUE)
})


test_that("smcfcs gives no warnings or errors with perfect prediction and brlogistic method
          for substantive model and covariate model", {
  expect_error({
    n <- 10
    x <- c(0,0,0,0,0,1,1,1,1,1)
    y <- c(0,0,0,0,0,1,1,1,1,1)
    x[c(5,10)] <- NA

    simData <- data.frame(x,y)
    rm(x,y)
    set.seed(62377)
    imps <- smcfcs(simData, smtype="brlogistic", smformula="y~x",
                   method=c("brlogreg", ""), m=1)
  }, NA)
})
