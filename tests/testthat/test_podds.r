library(smcfcs)
context("Testing imputation with podds using polr for ordered factors")

test_that("podds using polr runs, one covariate in polr model", {
  expect_error(
    {
      set.seed(1234)
      n <- 100

      z <- rnorm(n)
      xtilde <- 2*z + rlogis(n)

      x <- rep("a",n)
      x[(xtilde>(-1.5)) & (xtilde<0.25)] <- "b"
      x[(xtilde>(0.25)) & (xtilde<0.75)] <- "c"
      x[(xtilde>(0.75)) ] <- "d"
      table(x)

      y <- z+(x=="b")-(x=="c")+2*(x=="d")+rnorm(n)

      x[1*runif(n)<(exp(z)/(1+exp(z)))] <- NA

      simData <- data.frame(x=as.ordered(x),z=z,y=y)

      imps <- smcfcs(originaldata=simData,
        smtype="lm",
        smformula = "y~z+I(x=='b')+I(x=='c')+I(x=='d')",
        method = c("podds","", ""),
      )
    },
    NA
  )

})

test_that("podds using polr runs, two covariates in polr model", {
  expect_error(
    {
      set.seed(1234)
      n <- 100

      z1 <- rnorm(n)
      z2 <- z1+rnorm(n)
      xtilde <- 2*z1-z2 + rlogis(n)

      x <- rep("a",n)
      x[(xtilde>(-1.5)) & (xtilde<0.25)] <- "b"
      x[(xtilde>(0.25)) & (xtilde<0.75)] <- "c"
      x[(xtilde>(0.75)) ] <- "d"
      table(x)

      y <- z1+z2+(x=="b")-(x=="c")+2*(x=="d")+rnorm(n)

      x[1*runif(n)<(exp(z1)/(1+exp(z1)))] <- NA

      simData <- data.frame(x=as.ordered(x),z1=z1,z2=z2,y=y)

      imps <- smcfcs(originaldata=simData,
                     smtype="lm",
                     smformula = "y~z1+z2+I(x=='b')+I(x=='c')+I(x=='d')",
                     method = c("podds","","","")
      )
    },
    NA
  )
})


test_that("podds using polr is unbiased", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(2347633)
      n <- 100000

      z <- rnorm(n)
      xtilde <- 2*z + rlogis(n)

      x <- rep("a",n)
      x[(xtilde>(-1.5)) & (xtilde<0.25)] <- "b"
      x[(xtilde>(0.25)) & (xtilde<0.75)] <- "c"
      x[(xtilde>(0.75)) ] <- "d"
      table(x)

      y <- z+(x=="b")-(x=="c")+2*(x=="d")+rnorm(n)

      x[1*runif(n)<(exp(z)/(1+exp(z)))] <- NA

      simData <- data.frame(x=as.ordered(x),z=z,y=y)

      imps <- smcfcs(originaldata=simData,
                     smtype="lm",
                     smformula = "y~z+I(x=='b')+I(x=='c')+I(x=='d')",
                     method = c("podds","", ""),
      )

      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, lm(y ~ z+I(x=='b')+I(x=='c')+I(x=='d')))
      MIcombineRes <- summary(MIcombine(models))
      # 95% CI includes true value
      (MIcombineRes$`(lower`[5] < 2) & (MIcombineRes$`upper)`[5]>2)
    },
    TRUE
  )

})


test_that("podds using polr with two covariates is unbiased", {
  skip_on_cran()
  expect_equal(
    {
      set.seed(1762348012)
      n <- 100000

      z1 <- rnorm(n)
      z2 <- z1+rnorm(n)
      xtilde <- 2*z1-z2 + rlogis(n)

      x <- rep("a",n)
      x[(xtilde>(-1.5)) & (xtilde<0.25)] <- "b"
      x[(xtilde>(0.25)) & (xtilde<0.75)] <- "c"
      x[(xtilde>(0.75)) ] <- "d"
      table(x)

      y <- z1+z2+(x=="b")-(x=="c")+2*(x=="d")+rnorm(n)

      x[1*runif(n)<(exp(z1)/(1+exp(z1)))] <- NA

      simData <- data.frame(x=as.ordered(x),z1=z1,z2=z2,y=y)

      imps <- smcfcs(originaldata=simData,
                     smtype="lm",
                     smformula = "y~z1+z2+I(x=='b')+I(x=='c')+I(x=='d')",
                     method = c("podds","","","")
      )

      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, lm(y ~ z1+z2+I(x=='b')+I(x=='c')+I(x=='d')))
      MIcombineRes <- summary(MIcombine(models))
      # 95% CI includes true value
      (MIcombineRes$`(lower`[6] < 2) & (MIcombineRes$`upper)`[6]>2)
    },
    TRUE
  )
})
