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
  expect_equal({
    n <- 10000
    x <- rnorm(n)
    psxb <- x
    pspr <- exp(psxb)/(1+exp(psxb))
    z <- 1*(runif(n)<pspr)
    y <- x + z + rnorm(n)
    data <- data.frame(x,z,y)
    data$x[runif(n)<0.5] <- NA
    imps <- mips(data, omtype="lm", omformula="y~x+z",
                 method=c("norm","",""),psformula="z~x",m=1)
    imp1 <- imps$impDatasets[[1]]
    outmod <- lm(y~z+x, data=imp1)
    as.logical(abs(outmod$coef[2]-1)<0.05)
  }, TRUE)
})
