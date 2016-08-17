library(smcfcs)
context("Testing mips")

test_that("mips fails if no psformula is supplied", {
  expect_error({
    n <- 100
    x <- rnorm(n)
    psxb <- x
    pspr <- exp(psxb)/(1+exp(psxb))
    z <- 1*(runif(n)<pspr)
    y <- x + z + rnorm(n)
    data <- data.frame(x,z,y)
    data$x[runif(n)<0.5] <- NA
    imps <- mips(data, omtype="lm", omformula="y~x+z",
                 method=c("norm","",""))
  })
})

test_that("mips gives approximately unbiased estimates when it should", {
  testthat::skip_on_cran()
  set.seed(542)
  expect_equal({
    #we perform a small simulation study to check unbiasedness
    #first with a continuous confounder
    nSim <- 100
    ests <- array(0, dim=c(nSim,3))
    n <- 100
    m <- 5
    impEsts <- array(0, dim=c(m,3))

    for (sim in 1:nSim) {
      print(sim)
      x <- rnorm(n)
      psxb <- x
      pspr <- exp(psxb)/(1+exp(psxb))
      z <- 1*(runif(n)<pspr)
      y <- x + z + rnorm(n)
      data <- data.frame(x,z,y)
      data$x[runif(n)<0.5] <- NA
      imps <- mips(data, omtype="lm", omformula="y~x+z",
                   method=c("norm","",""),psformula="z~x",m=m)
      for (i in 1:m) {
        omFit <- lm(y~x+z, data=imps[[i]])
        impEsts[i,] <- coef(omFit)
      }
      ests[sim,] <- colMeans(impEsts)
    }
    as.logical(((mean(ests[,1])-0) <0.05) & ((mean(ests[,2])-1) <0.05) & ((mean(ests[,1])-1) <0.05))
  }, TRUE)

  expect_equal({
    #we perform a small simulation study to check unbiasedness
    #now with a binary confounder
    nSim <- 100
    ests <- array(0, dim=c(nSim,3))
    n <- 100
    m <- 5
    impEsts <- array(0, dim=c(m,3))

    for (sim in 1:nSim) {
      print(sim)
      x <- 1*(runif(n)<0.5)
      psxb <- x
      pspr <- exp(psxb)/(1+exp(psxb))
      z <- 1*(runif(n)<pspr)
      y <- x + z + rnorm(n)
      data <- data.frame(x,z,y)
      data$x[runif(n)<0.5] <- NA
      imps <- mips(data, omtype="lm", omformula="y~x+z",
                   method=c("logreg","",""),psformula="z~x",m=m)
      for (i in 1:m) {
        omFit <- lm(y~x+z, data=imps[[i]])
        impEsts[i,] <- coef(omFit)
      }
      ests[sim,] <- colMeans(impEsts)
    }
    as.logical(((mean(ests[,1])-0) <0.05) & ((mean(ests[,2])-1) <0.05) & ((mean(ests[,1])-1) <0.05))
  }, TRUE)
})
