library(smcfcs)
library(survival)
library(mitools)
context("Covariate measurement error testing")

test_that("Linear regression with cov. measurement error runs", {
  expect_error(
    {
      set.seed(1234)

      n <- 100 # sample size
      x <- rnorm(n)
      w1 <- x + rnorm(n)
      w2 <- x + rnorm(n)
      y <- x + rnorm(n)
      x <- rep(NA, n)

      simData <- data.frame(x, w1, w2, y)

      errMat <- matrix(0, nrow = 4, ncol = 4)
      errMat[1, c(2, 3)] <- 1

      imps <- smcfcs(simData,
        smtype = "lm", smformula = "y~x",
        method = c("latnorm", "", "", ""),
        errorProneMatrix = errMat, numit = 100, m = 1
      )
    },
    NA
  )
})

test_that("Linear regression with cov. measurement error is consistent", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(1234)

      n <- 10000 # sample size
      x <- rnorm(n)
      w1 <- x + rnorm(n)
      w2 <- x + rnorm(n)
      y <- x + rnorm(n)
      x <- rep(NA, n)

      simData <- data.frame(x, w1, w2, y)

      errMat <- matrix(0, nrow = 4, ncol = 4)
      errMat[1, c(2, 3)] <- 1

      imps <- smcfcs(simData,
        smtype = "lm", smformula = "y~x",
        method = c("latnorm", "", "", ""),
        errorProneMatrix = errMat, numit = 100, m = 1
      )
      as.logical(abs(coef(lm(y ~ x, data = imps$impDatasets[[1]]))[2] - 1) < 0.05)
    },
    TRUE
  )
})

test_that("Linear regression with cov. measurement error coverage is ok", {
  skip_on_cran()
  expect_equal(
    {
      library(mitools)
      set.seed(1234)
      nSim <- 1000
      n <- 1000
      ests <- array(0, dim = c(nSim, 2))

      for (i in 1:nSim) {
        x <- rnorm(n)
        w1 <- x + rnorm(n)
        w2 <- x + rnorm(n)
        w2[round(n * 0.1):n] <- NA
        y <- x + rnorm(n)
        x <- rep(NA, n)

        simData <- data.frame(x, w1, w2, y)

        errMat <- matrix(0, nrow = 4, ncol = 4)
        errMat[1, c(2, 3)] <- 1

        imps <- smcfcs(
          originaldata = simData,
          smtype = "lm",
          smformula = "y~x",
          method = c("latnorm", "", "", ""),
          errorProneMatrix = errMat, numit = 100, m = 5
        )
        impobj <- imputationList(imps$impDatasets)
        modelSummary <- summary(MIcombine(with(impobj, lm(y ~ x))))
        ests[i, ] <- c(modelSummary[2, 3], modelSummary[2, 4])
      }
      # check coverage is close to 95%
      abs(mean((ests[, 1] < 1) & (ests[, 2] > 1)) - 0.95) < (qnorm(0.999) * ((0.95 * 0.05) / nSim)^0.5)
    },
    TRUE
  )
})
