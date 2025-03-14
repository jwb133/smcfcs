library(smcfcs)
context("Testing imputation with restrictions for handling coarsening")

test_that("Restrictions runs and respects restrictions with coarsened factor covariate", {
  expect_equal(
    {
      set.seed(1234)
      n <- 100
      z <- rnorm(n)

      p1 <- 1/(1+exp(0+z)+exp(1+z))
      p2 <- exp(0+z)/(1+exp(0+z)+exp(1+z))
      p3 <- 1-p1-p2
      u <- runif(n)
      x <- rep("a",n)
      x[(u>p1) & (u<(p1+p2))] <- "b"
      x[(u>(p1+p2))] <- "c"

      y <- z+(x=="b")+2*(x=="c")+rnorm(n)

      # generate coarsening variable, where 0 indicates no coarsening
      # 1 indicates missing
      # 2 indicates coarsened variable=a if x=a or b/c if x=b or x=c
      # 3 indicates coarsened variable=b if x=b or a/c if x=a or x=c
      c <- sample(c(0, 1, 2, 3), size = n, replace = TRUE, prob = rep(0.25,4))
      xobs <- x
      xobs[c==1] <- "NA"
      xobs[(c==2) & ((x=="b") | (x=="c"))] <- "b/c"
      xobs[(c==3) & ((x=="a") | (x=="c"))] <- "a/c"
      table(xobs, useNA = "always")

      # set x to NA when not exactly observed
      x[c==1] <- NA
      x[(c==2) & ((x=="b") | (x=="c"))] <- NA
      x[(c==3) & ((x=="a") | (x=="c"))] <- NA

      simData <- data.frame(x=factor(x),xobs=xobs,z=z,y=y)

      restrictionsX = c("xobs = a/c ~ a + c",
                        "xobs = b/c ~ b + c")
      restrictions = append(list(restrictionsX), as.list(c("", "", "", "")))

      imps <- smcfcs(originaldata=simData,
        smtype="lm",
        smformula = "y~z+x",
        method = c("mlogit","", "", ""),
        restrictions = restrictions
      )

      # check imputations respect restrictions
      fail <- 0
      for (i in 1:length(imps$impDatasets)) {
        rowsToCheck <- subset(imps$impDatasets[[i]], xobs %in% c("a/c", "b/c"))
        rowsToCheck$valid <- with(rowsToCheck,
                               ifelse(xobs == "a/c", x %in% c("a", "c"),
                                      ifelse(xobs == "b/c", x %in% c("b", "c"), NA)))
        fail <- fail + sum(!rowsToCheck$valid)
      }
      fail
    },
    0
  )

})


test_that("Restrictions is unbiased with coarsened factor covariate", {
  skip_on_cran()
  expect_equal({
      set.seed(79968)
      n <- 100000
      z <- rnorm(n)

      p1 <- 1/(1+exp(0+z)+exp(1+z))
      p2 <- exp(0+z)/(1+exp(0+z)+exp(1+z))
      p3 <- 1-p1-p2
      u <- runif(n)
      x <- rep("a",n)
      x[(u>p1) & (u<(p1+p2))] <- "b"
      x[(u>(p1+p2))] <- "c"

      y <- z+(x=="b")+2*(x=="c")+rnorm(n)

      # generate coarsening variable, where 0 indicates no coarsening
      # 1 indicates missing
      # 2 indicates coarsened variable=a if x=a or b/c if x=b or x=c
      # 3 indicates coarsened variable=b if x=b or a/c if x=a or x=c
      c <- sample(c(0, 1, 2, 3), size = n, replace = TRUE, prob = rep(0.25,4))
      xobs <- x
      xobs[c==1] <- "NA"
      xobs[(c==2) & ((x=="b") | (x=="c"))] <- "b/c"
      xobs[(c==3) & ((x=="a") | (x=="c"))] <- "a/c"
      table(xobs, useNA = "always")

      # set x to NA when not exactly observed
      x[c==1] <- NA
      x[(c==2) & ((x=="b") | (x=="c"))] <- NA
      x[(c==3) & ((x=="a") | (x=="c"))] <- NA

      simData <- data.frame(x=factor(x),xobs=xobs,z=z,y=y)

      restrictionsX = c("xobs = a/c ~ a + c",
                        "xobs = b/c ~ b + c")
      restrictions = append(list(restrictionsX), as.list(c("", "", "", "")))

      imps <- smcfcs(originaldata=simData,
                     smtype="lm",
                     smformula = "y~z+x",
                     method = c("mlogit","", "", ""),
                     restrictions = restrictions
      )

      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, lm(y ~ z+x))
      MIcombineRes <- summary(MIcombine(models))
      # 95% CI includes true value
      (MIcombineRes$`(lower`[4] < 2) & (MIcombineRes$`upper)`[4]>2)
    },
    TRUE
  )

})


test_that("Restrictions runs respecting restrictions with coarsened ordered factor covariate", {
  expect_equal(
    {
      set.seed(1234)
      n <- 100

      z <- rnorm(n)
      xtilde <- 2*z + rlogis(n)

      x <- rep("a",n)
      x[(xtilde>(-1)) & (xtilde<1)] <- "b"
      x[(xtilde>(1)) ] <- "c"
      table(x)

      y <- z+(x=="b")+2*(x=="c")+rnorm(n)

      # generate coarsening variable, where 0 indicates no coarsening
      # 1 indicates missing
      # 2 indicates coarsened variable=a if x=a or b/c if x=b or x=c
      # 3 indicates coarsened variable=b if x=b or a/c if x=a or x=c
      c <- sample(c(0, 1, 2, 3), size = n, replace = TRUE, prob = rep(0.25,4))
      xobs <- x
      xobs[c==1] <- "NA"
      xobs[(c==2) & ((x=="b") | (x=="c"))] <- "b/c"
      xobs[(c==3) & ((x=="a") | (x=="c"))] <- "a/c"
      table(xobs, useNA = "always")

      # set x to NA when not exactly observed
      x[c==1] <- NA
      x[(c==2) & ((x=="b") | (x=="c"))] <- NA
      x[(c==3) & ((x=="a") | (x=="c"))] <- NA

      simData <- data.frame(x=as.ordered(x),xobs=xobs,z=z,y=y)

      restrictionsX = c("xobs = a/c ~ a + c",
                        "xobs = b/c ~ b + c")
      restrictions = append(list(restrictionsX), as.list(c("", "", "", "")))

      imps <- smcfcs(originaldata=simData,
                     smtype="lm",
                     smformula = "y~z+x",
                     method = c("podds","", "", ""),
                     restrictions = restrictions
      )

      # check imputations respect restrictions
      fail <- 0
      for (i in 1:length(imps$impDatasets)) {
        rowsToCheck <- subset(imps$impDatasets[[i]], xobs %in% c("a/c", "b/c"))
        rowsToCheck$valid <- with(rowsToCheck,
                                  ifelse(xobs == "a/c", x %in% c("a", "c"),
                                         ifelse(xobs == "b/c", x %in% c("b", "c"), NA)))
        fail <- fail + sum(!rowsToCheck$valid)
      }
      fail
    },
    0
  )

})



test_that("Restrictions is unbiased with coarsened ordered factor covariate", {
  expect_equal(
    {
      set.seed(79968)
      n <- 100000

      z <- rnorm(n)
      xtilde <- 2*z + rlogis(n)

      x <- rep("a",n)
      x[(xtilde>(-1)) & (xtilde<1)] <- "b"
      x[(xtilde>(1)) ] <- "c"
      table(x)

      y <- z+(x=="b")+2*(x=="c")+rnorm(n)

      # generate coarsening variable, where 0 indicates no coarsening
      # 1 indicates missing
      # 2 indicates coarsened variable=a if x=a or b/c if x=b or x=c
      # 3 indicates coarsened variable=b if x=b or a/c if x=a or x=c
      c <- sample(c(0, 1, 2, 3), size = n, replace = TRUE, prob = rep(0.25,4))
      xobs <- x
      xobs[c==1] <- "NA"
      xobs[(c==2) & ((x=="b") | (x=="c"))] <- "b/c"
      xobs[(c==3) & ((x=="a") | (x=="c"))] <- "a/c"
      table(xobs, useNA = "always")

      # set x to NA when not exactly observed
      x[c==1] <- NA
      x[(c==2) & ((x=="b") | (x=="c"))] <- NA
      x[(c==3) & ((x=="a") | (x=="c"))] <- NA

      simData <- data.frame(x=as.ordered(x),xobs=xobs,z=z,y=y)

      restrictionsX = c("xobs = a/c ~ a + c",
                        "xobs = b/c ~ b + c")
      restrictions = append(list(restrictionsX), as.list(c("", "", "", "")))

      imps <- smcfcs(originaldata=simData,
                     smtype="lm",
                     smformula = "y~z+x",
                     method = c("podds","", "", ""),
                     restrictions = restrictions
      )

      library(mitools)
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, lm(y ~ z+I(x=="b")+I(x=="c")))
      MIcombineRes <- summary(MIcombine(models))
      # 95% CI includes true value
      (MIcombineRes$`(lower`[4] < 2) & (MIcombineRes$`upper)`[4]>2)
    },
    TRUE
  )

})

