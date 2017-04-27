library(smcfcs)
library(survival)
library(mitools)
context("Case cohort testing")

test_that("Case cohort imputation runs and is approximately unbiased, binary covariate missing", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)

    #run simulation study
    nSims <- 100
    xEsts <- array(0, dim=nSims)
    n <- 10000

    for (sim in 1:nSims) {
      print(sim)
      z <- rnorm(n)
      x <- 1*(runif(n)<(exp(z)/(1+exp(z))))
      t <- -log(runif(n))/(0.01*exp(x+z))
      d <- 1*(t<2)
      t[d==0] <- 2
      x[(runif(n)<0.5)] <- NA

      fullcohortdata <- data.frame(t,d,x,z)
      fullcohortdata$in.subco <- 0
      #we sample a 10% subcohort
      fullcohortdata$in.subco[sample(n, size=n*0.1)] <- 1
      fullcohortdata$id <- 1:n

      ccdata <- fullcohortdata[(fullcohortdata$in.subco==1) | (fullcohortdata$d==1),]
      ccdata$entertime <- 0
      ccdata$entertime[ccdata$in.subco==0] <- ccdata$t[ccdata$in.subco==0] - 0.000001
      imps <- smcfcs.casecohort(ccdata, smformula="Surv(entertime, t, d)~x+z", sampfrac=0.1,
                                in.subco="in.subco", method=c("", "", "logreg", "", "", "", ""))
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, coxph(Surv(entertime,t,d)~x+z+cluster(id)))
      xEsts[sim] <- summary(MIcombine(models))[1,1]
    }
    print(mean(xEsts))
    abs(mean(xEsts-1))<0.1
  }, TRUE)
})

test_that("Case cohort imputation runs and is approximately unbiased, continuous covariate missing", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)

    #run simulation study
    nSims <- 100
    xEsts <- array(0, dim=nSims)
    n <- 10000

    for (sim in 1:nSims) {
      print(sim)
      z <- rnorm(n)
      x <- z+rnorm(n)
      t <- -log(runif(n))/(0.01*exp(x+z))
      d <- 1*(t<1)
      t[d==0] <- 1
      x[(runif(n)<0.5)] <- NA

      fullcohortdata <- data.frame(t,d,x,z)
      fullcohortdata$in.subco <- 0
      #we sample a 10% subcohort
      fullcohortdata$in.subco[sample(n, size=n*0.1)] <- 1
      fullcohortdata$id <- 1:n

      ccdata <- fullcohortdata[(fullcohortdata$in.subco==1) | (fullcohortdata$d==1),]
      ccdata$entertime <- 0
      ccdata$entertime[ccdata$in.subco==0] <- ccdata$t[ccdata$in.subco==0] - 0.000001
      imps <- smcfcs.casecohort(ccdata, smformula="Surv(entertime, t, d)~x+z", sampfrac=0.1,
                                in.subco="in.subco", method=c("", "", "norm", "", "", "", ""))
      impobj <- imputationList(imps$impDatasets)
      models <- with(impobj, coxph(Surv(entertime,t,d)~x+z+cluster(id)))
      xEsts[sim] <- summary(MIcombine(models))[1,1]
    }
    print(mean(xEsts))
    abs(mean(xEsts-1))<0.1
  }, TRUE)
})

