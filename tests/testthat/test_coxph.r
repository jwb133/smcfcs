library(smcfcs)
library(survival)
context("Cox proportional hazards model testing")

test_that("Cox imputation is approximately unbiased", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 10000
    z <- rnorm(n)
    x <- z+rnorm(n)
    t <- -log(runif(n))/(0.01*exp(x+z))
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x,z)

    imps <- smcfcs(simData, smtype="coxph", smformula="Surv(t, d)~x+z",
                   method=c("", "", "norm", ""))
    library(mitools)
    impobj <- imputationList(imps$impDatasets)
    models <- with(impobj, coxph(Surv(t,d)~x+z))
    abs(summary(MIcombine(models))[1,1]-1)<0.1
  }, TRUE)
})

test_that("Cox imputation works with only one covariate", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 10000
    x <- rnorm(n)
    t <- -log(runif(n))/(0.01*exp(x))
    d <- 1*(t<10)
    t[d==0] <- 10
    x[(runif(n)<0.5)] <- NA

    simData <- data.frame(t,d,x)

    imps <- smcfcs(simData, smtype="coxph", smformula="Surv(t, d)~x",
                   method=c("", "", "norm"))
    library(mitools)
    impobj <- imputationList(imps$impDatasets)
    models <- with(impobj, coxph(Surv(t,d)~x))
    abs(summary(MIcombine(models))[1,1]-1)<0.1
  }, TRUE)
})

test_that("Cox imputation with cov. measurement error", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)

    nSim <- 5
    ests <- array(0, dim=c(nSim,3))

    n<- 1000 #sample size
    n1<- 100 #size of calibration sample
    xzcorr <- 0.25
    kappa<-2 #weibull shape parameter

    reliability <- 0.5
    b1 <- 0.5
    lambda<-0.00164 #weibull scale parameter
    library(MASS)
    cov <- mvrnorm(n, mu=c(0,0), Sigma=array(c(1,xzcorr,xzcorr,1),dim=c(2,2)))

    x <- cov[,1]
    z <- cov[,2]

    u<-runif(n,0,1)
    t.event<-((-log(u)/lambda)*(kappa^(lambda+1))*exp(-(b1*x+b1*z)))^(1/kappa)
    d<-ifelse(t.event<10,1,0)
    t<-ifelse(d==1,t.event,10)
    mean(d)

    sigma_u_sq <- 1/reliability - 1
    w1 <- x+rnorm(n, sd=sigma_u_sq^0.5)
    w2 <- x+rnorm(n, sd=sigma_u_sq^0.5)
    w2[(n1+1):n] <- NA

    simData <- data.frame(t,d,x=NA,z,w1,w2)

    errMat <- matrix(0, nrow=6, ncol=6)
    errMat[3,c(5,6)] <- 1

    imps <- smcfcs(simData, smtype="coxph", smformula="Surv(t, d)~x+z",
                   method=c("", "", "latnorm","", "", ""),
                   errorProneMatrix=errMat,numit=100)
    library(mitools)
    impobj <- imputationList(imps$impDatasets)
    models <- with(impobj, coxph(Surv(t,d)~x+z))
    }
  }, TRUE)
})
