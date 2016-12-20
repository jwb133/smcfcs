library(smcfcs)
context("Case cohort testing")

test_that("Case cohort imputation runs", {
  expect_equal({
    skip_on_cran()
    n <- 1000
    x <- rnorm(n)
    t <- -log(runif(n))/(0.01*exp(x))
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.2)] <- NA

    fullcohortdata <- data.frame(t,d,x)
    fullcohortdata$in.subco <- 0
    fullcohortdata$in.subco[sample(n, size=n*0.1)] <- 1
    fullcohortdata$id <- 1:n

    ccdata <- fullcohortdata[(fullcohortdata$in.subco==1) | (fullcohortdata$d==1),]
    ccdata$entertime <- 0
    ccdata$entertime[ccdata$in.subco==0] <- ccdata$t[ccdata$in.subco==0] - 0.001
    imps <- smcfcs.casecohort(ccdata, smformula="Surv(entertime, t, d)~x", sampfrac=0.1,
                              in.subco="in.subco", method=c("", "", "norm", "", "", ""))
  }, 1)
})
