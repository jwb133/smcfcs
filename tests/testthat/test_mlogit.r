library(smcfcs)
context("Testing imputation with mlogit")

test_that("mlogit runs without error", {
  expect_no_error(
    {
      set.seed(1234)
      n <- 1000

      z <- rnorm(n)
      xb2 <- exp(-1+1*z)
      xb3 <- exp(0.5-0.75*z)
      pr2 <- xb2/(1+xb2+xb3)
      pr3 <- xb3/(1+xb2+xb3)
      pr1 <- 1-pr2-pr3
      u <- runif(n)
      x <- rep("a",n)
      x[(u>pr1) & (u>(pr1+pr2))] <- "c"
      x[(u>pr1) & (u<(pr1+pr2))] <- "b"

      y <- z+(x=="b")-(x=="c")+rnorm(n)

      x[1*runif(n)<(exp(z)/(1+exp(z)))] <- NA

      simData <- data.frame(x=as.factor(x),z=z,y=y)

      imps <- smcfcs(originaldata=simData,
        smtype="lm",
        smformula = "y~z+I(x=='b')+I(x=='c')",
        method = c("mlogit","", "")
      )
    },
  )
})

test_that("mlogit covariate imp is unbiased", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(2347633)
      n <- 100000

      z <- rnorm(n)
      xb2 <- exp(-1+1*z)
      xb3 <- exp(0.5-0.75*z)
      pr2 <- xb2/(1+xb2+xb3)
      pr3 <- xb3/(1+xb2+xb3)
      pr1 <- 1-pr2-pr3
      u <- runif(n)
      x <- rep("a",n)
      x[(u>pr1) & (u>(pr1+pr2))] <- "c"
      x[(u>pr1) & (u<(pr1+pr2))] <- "b"

      y <- z+(x=="b")-(x=="c")+rnorm(n)

      x[1*runif(n)<(exp(z)/(1+exp(z)))] <- NA

      simData <- data.frame(x=as.factor(x),z=z,y=y)

      imps <- smcfcs(originaldata=simData,
                     smtype="lm",
                     smformula = "y~z+I(x=='b')+I(x=='c')",
                     method = c("mlogit","", "")
      )

      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, lm(y ~ z+I(x=='b')+I(x=='c')))
      MIcombineRes <- summary(MIcombine(models))
      # 95% CI includes true value
      (MIcombineRes$`(lower`[4] < (-1)) & (MIcombineRes$`upper)`[4]>(-1))
    },
    TRUE
  )

})
