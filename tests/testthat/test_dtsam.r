library(smcfcs)
library(survival)
context("Discrete time survival analysis testing")

test_that("Full data dtsam is unbiased", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 100000
    x1 <- 1*(runif(n)<0.5)
    x2 <- x1+rnorm(n)
    T <- 10
    #define vector of intercepts
    alpha <- seq(-1,-0.1,0.1)
    beta <- c(1,-1)
    yMat <- array(0, dim=c(n,T))
    for (i in 1:T) {
      yMat[,i] <- 1*(runif(n)<expit(alpha[i]+beta[1]*x1+beta[2]*x2))
    }
    failtime <- apply(yMat,1,function(x) which(x==1)[1])
    #event indicator
    d <- rep(1,n)
    d[is.na(failtime)] <- 0
    failtime[is.na(failtime)] <- 10
    mean(d)

    simData <- data.frame(x1=x1, x2=x2, failtime=failtime, d=d)

    #fit dtsam model
    longData <- survSplit(Surv(failtime,d)~x1+x2, data=simData, cut=1:T)
    mod <- glm(d~-1+factor(tstart)+x1+x2, family="binomial", data=longData)
    summary(mod)
    ciLower <- coef(mod)-qnorm(0.999)*sqrt(diag(vcov(mod)))
    ciUpper <- coef(mod)+qnorm(0.999)*sqrt(diag(vcov(mod)))
    ciIncluded <- (ciLower<c(alpha,beta))*(ciUpper>c(alpha,beta))

    mean(tail(ciIncluded,length(beta)))==1
  }, TRUE)
})

test_that("MAR is unbiased, binary covariate missing", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 100000
    x1 <- 1*(runif(n)<0.5)
    x2 <- x1+rnorm(n)
    T <- 10
    #define vector of intercepts
    alpha <- seq(-1,-0.1,0.1)
    beta <- c(1,-1)
    yMat <- array(0, dim=c(n,T))
    for (i in 1:T) {
      yMat[,i] <- 1*(runif(n)<expit(alpha[i]+beta[1]*x1+beta[2]*x2))
    }
    failtime <- apply(yMat,1,function(x) which(x==1)[1])
    #event indicator
    d <- rep(1,n)
    d[is.na(failtime)] <- 0
    failtime[is.na(failtime)] <- 10
    mean(d)
    #add in some additional completely random censoring
    cTime <- ceiling(10*runif(n))
    t <- failtime
    t[cTime<failtime] <- cTime[cTime<failtime]
    d[cTime<failtime] <- 0
    mean(d)

    simData <- data.frame(x1=x1, x2=x2, failtime=failtime, d=d)
    simData$x1[runif(n)<expit((x2-mean(x2)/sd(x2)))] <- NA
    #simData$x2[runif(n)<0.25] <- NA

    #impute using smcfcs.dtsam
    M <- 5
    imps <- smcfcs.dtsam(simData, "Surv(failtime,d)~x1+x2",
                   method=c("logreg","", "", ""),m=M)

    #fit dtsam model to each dataset
    ests <- vector(mode = "list", length = M)
    vars <- vector(mode = "list", length = M)
    for (i in 1:M) {
      longData <- survSplit(Surv(failtime,d)~x1+x2, data=imps$impDatasets[[i]],
                            cut=1:T)
      mod <- glm(d~-1+factor(tstart)+x1+x2, family="binomial", data=longData)
      ests[[i]] <- coef(mod)
      vars[[i]] <- diag(vcov(mod))
    }

    rubin <- mitools::MIcombine(ests,vars)

    ciLower <- rubin$coefficients-qt(0.999,df=rubin$df)*sqrt(diag(rubin$variance))
    ciUpper <- rubin$coefficients+qt(0.999,df=rubin$df)*sqrt(diag(rubin$variance))
    ciIncluded <- (ciLower<c(alpha,beta))*(ciUpper>c(alpha,beta))

    mean(tail(ciIncluded,length(beta)))==1
  }, TRUE)
})

test_that("MAR is unbiased, cts covariate missing", {
  skip_on_cran()
  expect_equal({
    set.seed(1234612)
    n <- 1000

    x1 <- 1*(runif(n)<0.5)
    x2 <- x1+rnorm(n)
    T <- 10
    #define vector of intercepts
    alpha <- seq(-1,-0.1,0.1)
    beta <- c(1,-1)
    yMat <- array(0, dim=c(n,T))
    for (i in 1:T) {
      yMat[,i] <- 1*(runif(n)<expit(alpha[i]+beta[1]*x1+beta[2]*x2))
    }
    failtime <- apply(yMat,1,function(x) which(x==1)[1])
    #event indicator
    d <- rep(1,n)
    d[is.na(failtime)] <- 0
    failtime[is.na(failtime)] <- 10
    mean(d)
    #add in some additional completely random censoring
    cTime <- ceiling(10*runif(n))
    t <- failtime
    t[cTime<failtime] <- cTime[cTime<failtime]
    d[cTime<failtime] <- 0
    mean(d)

    simData <- data.frame(x1=x1, x2=x2, failtime=failtime, d=d)
    #simData$x2[runif(n)<expit((x1-mean(x1)/sd(x1)))] <- NA
    simData$x2[runif(n)<0.5] <- NA

    #impute using smcfcs.dtsam
    M <- 5
    imps <- smcfcs.dtsam(simData, "Surv(failtime,d)~x1+x2",
                         method=c("","norm", "", ""),m=M)

    #fit dtsam model to each dataset
    ests <- vector(mode = "list", length = M)
    vars <- vector(mode = "list", length = M)
    for (i in 1:M) {
      longData <- survSplit(Surv(failtime,d)~x1+x2, data=imps$impDatasets[[i]],
                            cut=1:T)
      mod <- glm(d~-1+factor(tstart)+x1+x2, family="binomial", data=longData)
      ests[[i]] <- coef(mod)
      vars[[i]] <- diag(vcov(mod))
    }

    rubin <- mitools::MIcombine(ests,vars)

    ciLower <- rubin$coefficients-qt(0.999,df=rubin$df)*sqrt(diag(rubin$variance))
    ciUpper <- rubin$coefficients+qt(0.999,df=rubin$df)*sqrt(diag(rubin$variance))
    ciIncluded <- (ciLower<c(alpha,beta))*(ciUpper>c(alpha,beta))

    mean(tail(ciIncluded,length(beta)))==1
  }, TRUE)
})

test_that("MCAR is unbiased, binary covariate missing, linear time effect", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 100000
    x1 <- 1*(runif(n)<0.5)
    x2 <- x1+rnorm(n)
    T <- 10
    #define vector of intercepts
    alpha <- c(-2,0.1)
    beta <- c(1,-1)
    yMat <- array(0, dim=c(n,T))
    for (i in 1:T) {
      yMat[,i] <- 1*(runif(n)<expit(alpha[1]+alpha[2]*(i-1)+beta[1]*x1+beta[2]*x2))
    }
    failtime <- apply(yMat,1,function(x) which(x==1)[1])
    #event indicator
    d <- rep(1,n)
    d[is.na(failtime)] <- 0
    failtime[is.na(failtime)] <- 10
    mean(d)

    simData <- data.frame(x1=x1, x2=x2, failtime=failtime, d=d)

    #fit dtsam model to full dataset
    longData <- survSplit(Surv(failtime,d)~x1+x2, data=simData, cut=1:T)
    mod <- glm(d~tstart+x1+x2, family="binomial", data=longData)
    summary(mod)

    simData$x1[runif(n)<0.5] <- NA
    #simData$x2[runif(n)<0.25] <- NA

    #impute using smcfcs.dtsam
    M <- 5
    imps <- smcfcs.dtsam(simData, "Surv(failtime,d)~x1+x2",timeEffects="linear",
                         method=c("logreg","", "", ""),m=M)

    #fit dtsam model to each dataset
    ests <- vector(mode = "list", length = M)
    vars <- vector(mode = "list", length = M)
    for (i in 1:M) {
      longData <- survSplit(Surv(failtime,d)~x1+x2, data=imps$impDatasets[[i]],
                            cut=1:T)
      mod <- glm(d~tstart+x1+x2, family="binomial", data=longData)
      ests[[i]] <- coef(mod)
      vars[[i]] <- diag(vcov(mod))
    }

    rubin <- mitools::MIcombine(ests,vars)

    ciLower <- rubin$coefficients-qt(0.999,df=rubin$df)*sqrt(diag(rubin$variance))
    ciUpper <- rubin$coefficients+qt(0.999,df=rubin$df)*sqrt(diag(rubin$variance))
    ciIncluded <- (ciLower<c(alpha,beta))*(ciUpper>c(alpha,beta))

    mean(ciIncluded)==1
  }, TRUE)
})

test_that("MCAR is unbiased, binary covariate missing, quadratic time effect", {
  skip_on_cran()
  expect_equal({
    set.seed(1234)
    n <- 100000
    x1 <- 1*(runif(n)<0.5)
    x2 <- x1+rnorm(n)
    T <- 10
    #define vector of intercepts
    alpha <- c(-3,0.1,0.02)
    beta <- c(1,-1)
    yMat <- array(0, dim=c(n,T))
    for (i in 1:T) {
      yMat[,i] <- 1*(runif(n)<expit(alpha[1]+alpha[2]*(i-1)+alpha[3]*(i-1)^2+beta[1]*x1+beta[2]*x2))
    }
    failtime <- apply(yMat,1,function(x) which(x==1)[1])
    #event indicator
    d <- rep(1,n)
    d[is.na(failtime)] <- 0
    failtime[is.na(failtime)] <- 10
    mean(d)

    simData <- data.frame(x1=x1, x2=x2, failtime=failtime, d=d)

    #fit dtsam model to full data
    longData <- survSplit(Surv(failtime,d)~x1+x2, data=simData, cut=1:T)
    mod <- glm(d~tstart+I(tstart^2)+x1+x2, family="binomial", data=longData)
    summary(mod)

    simData$x1[runif(n)<0.5] <- NA
    #simData$x2[runif(n)<0.25] <- NA

    #impute using smcfcs.dtsam
    M <- 5
    imps <- smcfcs.dtsam(simData, "Surv(failtime,d)~x1+x2",timeEffects="quad",
                         method=c("logreg","", "", ""),m=M)

    #fit dtsam model to each dataset
    ests <- vector(mode = "list", length = M)
    vars <- vector(mode = "list", length = M)
    for (i in 1:M) {
      longData <- survSplit(Surv(failtime,d)~x1+x2, data=imps$impDatasets[[i]],
                            cut=1:T)
      mod <- glm(d~tstart+I(tstart^2)+x1+x2, family="binomial", data=longData)
      ests[[i]] <- coef(mod)
      vars[[i]] <- diag(vcov(mod))
    }

    rubin <- mitools::MIcombine(ests,vars)

    ciLower <- rubin$coefficients-qt(0.999,df=rubin$df)*sqrt(diag(rubin$variance))
    ciUpper <- rubin$coefficients+qt(0.999,df=rubin$df)*sqrt(diag(rubin$variance))
    ciIncluded <- (ciLower<c(alpha,beta))*(ciUpper>c(alpha,beta))

    mean(ciIncluded)==1
  }, TRUE)
})
