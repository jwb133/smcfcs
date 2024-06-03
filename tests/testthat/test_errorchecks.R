library(smcfcs)
library(survival)
context("Error trap testing")

test_that("Checking outcome model check for logistic models", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- as.factor(1 * (runif(n) < 0.5))
    x[(runif(n) < 0.5)] <- NA

    simData <- data.frame(x, y)

    imps <- smcfcs(simData,
      smtype = "logistic", smformula = "y~x",
      method = c("norm", "")
    )
  })
})

test_that("Checking outcome model check for logistic models", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- 1 + 1 * (runif(n) < 0.5)
    x[(runif(n) < 0.5)] <- NA

    simData <- data.frame(x, y)

    imps <- smcfcs(simData,
      smtype = "logistic", smformula = "y~x",
      method = c("norm", "")
    )
  })
})

test_that("Checking outcome model check for logistic models", {
  expect_output({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- 1 * (runif(n) < 0.5)
    x[(runif(n) < 0.5)] <- NA

    simData <- data.frame(x, y)

    imps <- smcfcs(simData,
      smtype = "logistic", smformula = "y~x",
      method = c("norm", "")
    )
  })
})

test_that("Checking error trap for method statement 1", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- y <- x1 + x2 + rnorm(n)
    x1[(runif(n) < 0.5)] <- NA

    simData <- data.frame(x1, x2, y)

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x1+x2",
      method = c("", "", "")
    )
  })
})

test_that("Checking error trap for method statement 2", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- y <- x1 + x2 + rnorm(n)
    x1[(runif(n) < 0.5)] <- NA

    simData <- data.frame(x1, x2, y)

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x1+x2",
      method = c("", "norm", "")
    )
  })
})

test_that("Checking measurement error error checks 1", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    w1 <- x + rnorm(n)
    w2 <- x + rnorm(n)
    y <- y <- x + rnorm(n)
    x <- rep(NA, n)

    simData <- data.frame(x, w1, w2, y)

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x",
      method = c("latnorm", "", "", "")
    )
  })
})

test_that("Checking measurement error error checks 2", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    w1 <- x + rnorm(n)
    y <- y <- x + rnorm(n)
    x <- rep(NA, n)

    simData <- data.frame(x, w1, y)
    errMat <- array(0, dim = c(3, 3))

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x",
      method = c("latnorm", "", ""), errorProneMatrix = errMat
    )
  })
})

test_that("Checking measurement error error checks 3", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    w1 <- x + rnorm(n)
    y <- y <- x + rnorm(n)
    x <- rep(NA, n)

    simData <- data.frame(x, w1, y)
    errMat <- array(0, dim = c(3, 3))
    errMat[1, 1] <- 1

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x",
      method = c("latnorm", "", ""), errorProneMatrix = errMat
    )
  })
})

test_that("Checking measurement error error checks 4", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    w1 <- x + rnorm(n)
    w2 <- x + rnorm(n)
    y <- y <- x + rnorm(n)
    x <- rep(NA, n)

    simData <- data.frame(x, w1, w2, y)
    errMat <- array(0, dim = c(4, 4))
    errMat[1, 2] <- 1
    errMat[1, 3] <- 1
    errMat[4, 2] <- 1

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x",
      method = c("latnorm", "", "", ""), errorProneMatrix = errMat
    )
  })
})

test_that("Checking measurement error error checks 5", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    w1 <- x + rnorm(n)
    w2 <- x + rnorm(n)
    y <- y <- x + rnorm(n)
    x <- rep(NA, n)

    simData <- data.frame(x, w1, w2, y)
    errMat <- array(0, dim = c(4, 4))
    errMat[1, 2] <- 1
    errMat[1, 3] <- 1
    errMat[4, 2] <- 2

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x",
      method = c("latnorm", "", "", ""), errorProneMatrix = errMat
    )
  })
})


test_that("Checking measurement error error checks 6", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    w1 <- x + rnorm(n)
    w2 <- x + rnorm(n)
    y <- y <- x + rnorm(n)
    x <- rep(NA, n)

    simData <- data.frame(x, w1, w2, y)
    errMat <- array(0, dim = c(5, 5))
    errMat[1, 2] <- 1
    errMat[1, 3] <- 1

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x",
      method = c("latnorm", "", "", ""), errorProneMatrix = errMat
    )
  })
})

test_that("Checking measurement error error checks 7", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    w1 <- x + rnorm(n)
    w2 <- x + rnorm(n)
    y <- y <- x + rnorm(n)
    x <- rep(NA, n)
    z <- rnorm(n)
    z[1:50] <- NA

    simData <- data.frame(w1, w2, y, z)
    errMat <- array(0, dim = c(4, 4))
    errMat[3, 1] <- 1
    errMat[3, 2] <- 1

    imps <- smcfcs(simData,
      smtype = "lm", smformula = "y~x",
      method = c("", "", "", "norm"), errorProneMatrix = errMat
    )
  })
})

test_that("Cox imputation fails if event indicator is not coded right", {
  expect_error({
    set.seed(1234)
    n <- 1000
    z <- rnorm(n)
    x <- z + rnorm(n)
    t <- -log(runif(n)) / (1 * exp(x + z))
    d <- 1 * (t < 10) + 1
    t[d == 1] <- 10
    x[(runif(n) < 0.5)] <- NA

    simData <- data.frame(t, d, x, z)

    imps <- smcfcs(simData,
      smtype = "coxph", smformula = "Surv(t, d)~x+z",
      method = c("", "", "norm", "")
    )
  })
})

test_that("Cox error check doesn't fail when events are 1 1 1 0 0 0", {
  expect_error(
    {
      set.seed(1234)
      n <- 100
      z <- rnorm(n)
      x <- z + rnorm(n)
      t <- runif(n)
      d <- c(rep(1, n / 2), rep(0, n / 2))
      x[(runif(n) < 0.5)] <- NA

      simData <- data.frame(t, d, x, z)

      imps <- smcfcs(simData,
        smtype = "coxph", smformula = "Surv(t, d)~x+z",
        method = c("", "", "norm", "")
      )
    },
    NA
  )
})
